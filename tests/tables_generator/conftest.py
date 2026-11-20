"""pytest plugin producing an HTML report of the table-generator tests.

After every run of ``pytest tests/tables_generator`` a self-contained HTML file is
written (default: tests/tables_generator/report/test_report.html, change it with
``--tg-report=<path>``, disable with ``--tg-report=``). For every test it shows:
  * what is tested (first paragraph of the test docstring),
  * the result (passed / failed / known bug / skipped),
  * the reason of a known bug (xfail) or the error message of a failure,
  * the figures and tables the test attached through the ``report`` fixture.

Figures need matplotlib; without it the report is written without figures.
"""
from __future__ import annotations

import base64
import contextlib
import datetime as dt
import html
import io
import platform
import shutil
import subprocess
from collections import OrderedDict
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
DEFAULT_REPORT = HERE / "report" / "test_report.html"

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except Exception:  # pragma: no cover - matplotlib is optional
    plt = None


# ---------------------------------------------------------------------------------------------
# Collected data
# ---------------------------------------------------------------------------------------------
class _Null:
    """Stand-in for matplotlib axes when matplotlib is missing: accepts and ignores every call."""

    def __getattr__(self, _):
        return lambda *a, **k: _Null()

    def __getitem__(self, _):
        return self

    def __iter__(self):
        return iter(())


class TestRecord:
    def __init__(self, item):
        self.nodeid = item.nodeid
        self.name = item.name
        self.func = getattr(item, "originalname", item.name)
        cls = getattr(item, "cls", None)
        self.section = (cls.__doc__ or cls.__name__).strip().splitlines()[0] if cls else "Other tests"
        self.section_doc = _paragraphs(cls.__doc__)[1:] if cls else []
        self.doc = _paragraphs(getattr(item.function, "__doc__", None))
        self.outcome = "not run"
        self.detail = ""
        self.duration = 0.0
        self.figures: list[tuple[str, str]] = []  # (caption, base64 png)
        self.tables: list[tuple[str, list, list]] = []
        self.notes: list[str] = []


def _paragraphs(doc):
    if not doc:
        return []
    import inspect
    return [" ".join(p.split()) for p in inspect.cleandoc(doc).split("\n\n") if p.strip()]


RECORDS: "OrderedDict[str, TestRecord]" = OrderedDict()


class Report:
    """Fixture object used by tests to attach figures, tables and notes to the HTML report."""

    def __init__(self, record: TestRecord):
        self.record = record

    @contextlib.contextmanager
    def figure(self, caption: str, nrows: int = 1, ncols: int = 1, width: float = 5.0, height: float = 3.6):
        """with report.figure("caption", ncols=2) as axes: axes[0].plot(...)"""
        if plt is None:
            yield _Null()
            return
        fig, axes = plt.subplots(nrows, ncols, figsize=(width * ncols, height * nrows), squeeze=False)
        flat = axes.ravel()
        try:
            yield flat if flat.size > 1 else flat[0]
        finally:
            for ax in flat:
                ax.grid(True, alpha=0.3)
                if ax.get_legend_handles_labels()[0]:
                    ax.legend(fontsize=8)
            fig.tight_layout()
            buf = io.BytesIO()
            fig.savefig(buf, format="png", dpi=90)
            plt.close(fig)
            self.record.figures.append((caption, base64.b64encode(buf.getvalue()).decode()))

    def table(self, caption: str, headers: list, rows: list):
        self.record.tables.append((caption, headers, rows))

    def note(self, text: str):
        self.record.notes.append(text)


# ---------------------------------------------------------------------------------------------
# Hooks
# ---------------------------------------------------------------------------------------------
def pytest_addoption(parser):
    parser.addoption("--tg-report", default=str(DEFAULT_REPORT),
                     help="HTML report of the table-generator tests (empty string to disable)")


def _record(item) -> TestRecord:
    if item.nodeid not in RECORDS:
        RECORDS[item.nodeid] = TestRecord(item)
    return RECORDS[item.nodeid]


def pytest_collection_modifyitems(session, config, items):
    for item in items:
        if Path(item.fspath).resolve().parent == HERE:
            _record(item)


@pytest.fixture
def report(request):
    return Report(_record(request.node))


@pytest.hookimpl(hookwrapper=True)
def pytest_runtest_makereport(item, call):
    outcome = yield
    rep = outcome.get_result()
    if item.nodeid not in RECORDS:
        return
    rec = RECORDS[item.nodeid]
    rec.duration += rep.duration
    if rep.when == "setup":
        if rep.skipped:
            rec.outcome, rec.detail = "skipped", _skip_reason(rep)
        elif rep.failed:
            rec.outcome, rec.detail = "error", _short(rep.longreprtext)
    elif rep.when == "call":
        if hasattr(rep, "wasxfail"):
            rec.outcome = "known bug" if rep.skipped else "bug fixed?"
            rec.detail = rep.wasxfail
            if rep.skipped and call.excinfo is not None:
                rec.notes.append("Observed failure: " + _short(str(call.excinfo.value), 6))
        elif rep.passed:
            rec.outcome = "passed"
        elif rep.failed:
            text = rep.longreprtext
            if "XPASS(strict)" in text:
                rec.outcome = "bug fixed?"
                rec.detail = ("This test was marked as a known bug but now passes: the bug looks fixed, "
                              "remove the xfail marker. " + text.split("XPASS(strict)")[-1].strip(" ]"))
            else:
                rec.outcome, rec.detail = "failed", _short(text)
        elif rep.skipped:
            rec.outcome, rec.detail = "skipped", _skip_reason(rep)


def _skip_reason(rep):
    lr = rep.longrepr
    return lr[2] if isinstance(lr, tuple) and len(lr) == 3 else str(lr)


def _short(text, n=15):
    lines = [l for l in str(text).splitlines() if l.strip()]
    errs = [l for l in lines if l.lstrip().startswith("E ")]
    return "\n".join((errs or lines)[-n:])


def pytest_terminal_summary(terminalreporter, exitstatus, config):
    path = config.getoption("--tg-report")
    if not path or not any(r.outcome != "not run" for r in RECORDS.values()):
        return
    path = Path(path).resolve()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(_render())
    terminalreporter.write_sep("-", "table-generator test report")
    terminalreporter.write_line(f"open {path}")


# ---------------------------------------------------------------------------------------------
# HTML rendering
# ---------------------------------------------------------------------------------------------
STATUS = {  # label, css class, explanation
    "passed": ("PASSED", "ok", "works as expected"),
    "failed": ("FAILED", "bad", "unexpected failure: something is broken"),
    "error": ("ERROR", "bad", "the test could not be set up"),
    "known bug": ("KNOWN BUG", "warn", "documented bug, not fixed yet (xfail)"),
    "bug fixed?": ("BUG FIXED?", "info", "a known-bug test now passes: remove its xfail marker"),
    "skipped": ("SKIPPED", "skip", "not run (missing tool)"),
    "not run": ("NOT RUN", "skip", "deselected"),
}


def _env():
    def run(cmd):
        try:
            return subprocess.run(cmd, capture_output=True, text=True, cwd=HERE, timeout=10).stdout.strip()
        except Exception:
            return "?"
    from harness import lapack_dir
    gf = run(["gfortran", "--version"]).splitlines()[0] if shutil.which("gfortran") else "not found"
    return [("Date", dt.datetime.now().strftime("%Y-%m-%d %H:%M")),
            ("Git", f"{run(['git', 'rev-parse', '--abbrev-ref', 'HEAD'])} @ {run(['git', 'rev-parse', '--short', 'HEAD'])}"
                    + (" (with uncommitted changes)" if run(["git", "status", "--porcelain", "--", "../.."]) else "")),
            ("gfortran", gf), ("LAPACK", str(lapack_dir() or "not found")),
            ("Python", platform.python_version()), ("matplotlib", "yes" if plt else "no (no figures)")]


def _badge(outcome, count=None):
    label, css, _ = STATUS[outcome]
    return f'<span class="badge {css}">{label}{f" &times;{count}" if count else ""}</span>'


def _render():
    e = html.escape
    counts = {}
    for r in RECORDS.values():
        counts[r.outcome] = counts.get(r.outcome, 0) + 1
    sections: "OrderedDict[str, OrderedDict[str, list[TestRecord]]]" = OrderedDict()
    for r in RECORDS.values():
        sections.setdefault(r.section, OrderedDict()).setdefault(r.func, []).append(r)

    out = [f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Table generator tests</title><style>{CSS}</style></head><body><main>
<h1>Table generator &mdash; test report</h1>
<p class="lead">Automatic checks of <code>tables_generator/</code> (Tmat: scattering per diameter,
TmatInt: integration over the size distribution) and of <code>generate_tables.sh</code>.
Each card says what is checked, the result and, for physics checks, a plot of what was compared.</p>
<div class="summary">"""]
    for k in STATUS:
        if counts.get(k):
            out.append(f'<div class="tile {STATUS[k][1]}"><b>{counts[k]}</b><span>{STATUS[k][0].lower()}</span>'
                       f'<small>{STATUS[k][2]}</small></div>')
    out.append('</div><table class="env">' + "".join(f"<tr><th>{e(k)}</th><td>{e(v)}</td></tr>"
                                                   for k, v in _env()) + "</table>")
    out.append('<nav><b>Sections</b><ol>' + "".join(f'<li><a href="#s{i}">{e(s)}</a></li>'
                                                     for i, s in enumerate(sections)) + "</ol></nav>")

    for i, (section, funcs) in enumerate(sections.items()):
        first = next(iter(funcs.values()))[0]
        out.append(f'<section id="s{i}"><h2>{e(section)}</h2>')
        out += [f"<p>{e(p)}</p>" for p in first.section_doc]
        for func, recs in funcs.items():
            out.append(_render_card(func, recs))
        out.append("</section>")
    out.append("</main></body></html>")
    return "\n".join(out)


def _render_card(func, recs):
    e = html.escape
    rec0 = recs[0]
    per = {}
    for r in recs:
        per[r.outcome] = per.get(r.outcome, 0) + 1
    worst = next((o for o in ("failed", "error", "bug fixed?", "known bug", "skipped", "passed", "not run")
                  if o in per), "not run")
    badges = " ".join(_badge(o, n if len(recs) > 1 else None) for o, n in per.items())
    title = func.replace("test_", "", 1).replace("_", " ")
    out = [f'<details class="card {STATUS[worst][1]}" {"open" if worst != "passed" or any(r.figures for r in recs) else ""}>'
           f'<summary><span class="title">{e(title)}</span>{badges}</summary>']
    if rec0.doc:
        out.append(f'<p class="what"><b>What is checked:</b> {e(rec0.doc[0])}</p>')
        out += [f"<p>{e(p)}</p>" for p in rec0.doc[1:]]
    for r in recs:
        label = r.name[len(func):].strip("[]") if len(recs) > 1 else ""
        show = len(recs) == 1 or r.outcome != "passed" or r.figures or r.tables
        if not show:
            continue
        if label:
            out.append(f'<h4>{_badge(r.outcome)} <code>{e(label)}</code></h4>')
        if r.detail:
            cls = "reason" if r.outcome in ("known bug", "skipped", "bug fixed?") else "error"
            head = {"known bug": "Known bug", "skipped": "Skipped", "bug fixed?": "Note"}.get(r.outcome, "Error")
            out.append(f'<div class="{cls}"><b>{head}:</b> <pre>{e(r.detail)}</pre></div>')
        for n in r.notes:
            out.append(f'<p class="note">{e(n)}</p>')
        for cap, headers, rows in r.tables:
            out.append(f'<table class="data"><caption>{e(cap)}</caption><tr>'
                       + "".join(f"<th>{e(str(h))}</th>" for h in headers) + "</tr>"
                       + "".join("<tr>" + "".join(f"<td>{e(str(c))}</td>" for c in row) + "</tr>" for row in rows)
                       + "</table>")
        for cap, b64 in r.figures:
            out.append(f'<figure><img alt="{e(cap)}" src="data:image/png;base64,{b64}">'
                       f'<figcaption>{e(cap)}</figcaption></figure>')
    if len(recs) > 1:
        npass = per.get("passed", 0)
        out.append(f'<p class="note">{npass} of {len(recs)} cases passed'
                   + ("; passing cases are not listed." if npass and npass < len(recs) else ".") + "</p>")
    out.append(f'<p class="meta"><code>{e(rec0.nodeid.split("[")[0])}</code> &middot; '
               f'{sum(r.duration for r in recs):.2f} s</p></details>')
    return "".join(out)


CSS = """
:root{--bg:#fafaf9;--fg:#1c1917;--muted:#57534e;--card:#fff;--line:#e7e5e4;
--ok:#15803d;--okbg:#dcfce7;--bad:#b91c1c;--badbg:#fee2e2;--warn:#a16207;--warnbg:#fef3c7;
--info:#1d4ed8;--infobg:#dbeafe;--skip:#57534e;--skipbg:#f5f5f4}
@media (prefers-color-scheme:dark){:root{--bg:#1c1917;--fg:#f5f5f4;--muted:#a8a29e;--card:#292524;--line:#44403c;
--okbg:#14532d;--ok:#bbf7d0;--badbg:#7f1d1d;--bad:#fecaca;--warnbg:#713f12;--warn:#fde68a;--infobg:#1e3a8a;
--info:#bfdbfe;--skipbg:#44403c;--skip:#e7e5e4}figure img{background:#fff}}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--fg);font:15px/1.5 system-ui,sans-serif}
main{max-width:1100px;margin:auto;padding:24px 16px}h1{margin:.2em 0}.lead{color:var(--muted)}
h2{margin-top:2em;border-bottom:2px solid var(--line);padding-bottom:.2em}
.summary{display:flex;flex-wrap:wrap;gap:10px;margin:16px 0}
.tile{flex:1 1 150px;border-radius:10px;padding:10px 14px;display:flex;flex-direction:column}
.tile b{font-size:28px}.tile small{opacity:.8}
.ok{--c:var(--ok);--b:var(--okbg)}.bad{--c:var(--bad);--b:var(--badbg)}.warn{--c:var(--warn);--b:var(--warnbg)}
.info{--c:var(--info);--b:var(--infobg)}.skip{--c:var(--skip);--b:var(--skipbg)}
.tile,.badge{background:var(--b);color:var(--c)}
.badge{font-size:12px;font-weight:600;border-radius:999px;padding:2px 9px;margin-left:6px;white-space:nowrap}
table.env{border-collapse:collapse;font-size:13px;margin:8px 0}table.env th{text-align:left;padding-right:12px;color:var(--muted)}
nav{background:var(--card);border:1px solid var(--line);border-radius:10px;padding:8px 16px;margin:16px 0}
.card{background:var(--card);border:1px solid var(--line);border-left:5px solid var(--c);border-radius:8px;
margin:10px 0;padding:6px 14px}
summary{cursor:pointer;padding:6px 0;display:flex;flex-wrap:wrap;align-items:center;gap:4px}
.title{font-weight:600;flex:1 1 300px}.what{margin-top:4px}.meta,.note{color:var(--muted);font-size:13px}
pre{white-space:pre-wrap;margin:4px 0;font-size:12.5px}.reason,.error{border-radius:6px;padding:6px 10px;margin:8px 0}
.reason{background:var(--warnbg)}.error{background:var(--badbg)}
figure{margin:10px 0}figure img{max-width:100%;height:auto;border:1px solid var(--line);border-radius:6px}
figcaption{font-size:13px;color:var(--muted)}
table.data{border-collapse:collapse;font-size:13px;margin:8px 0;display:block;overflow-x:auto}
table.data td,table.data th{border:1px solid var(--line);padding:3px 8px}caption{text-align:left;color:var(--muted)}
"""
