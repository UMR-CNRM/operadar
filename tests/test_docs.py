#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""The documentation (docs/wiki/, published as the GitHub wiki, and the README files) must stay consistent with
the code: configuration options, operadar() arguments, example configurations, files and scripts named in the
pages, links between wiki pages and supported Python versions. When one of these tests fails after a code
change, update the page named in the message."""
import inspect
import re
import tomllib
from pathlib import Path

import pytest

from operadar.forward_operator import Config, operadar

ROOT = Path(__file__).resolve().parents[1]
WIKI = ROOT / "docs" / "wiki"
PAGES = sorted(WIKI.glob("*.md"))
DOCS = PAGES + [ROOT / "README.md", ROOT / "plot_tools" / "README.md"]
WIKI_URL = "https://github.com/UMR-CNRM/operadar/wiki"
BLOB_URL = "https://github.com/UMR-CNRM/operadar/blob/master/"


def read(path):
    return path.read_text(encoding="utf-8")


def first_column_names(page, header):
    """Names in backquotes of the first column of the markdown table that follows the line `header`."""
    lines = read(WIKI / page).splitlines()
    start = next(i for i, line in enumerate(lines) if line.strip() == header)
    names, in_table = [], False
    for line in lines[start + 1:]:
        if line.startswith("|"):
            in_table = True
            m = re.match(r"\|\s*`([^`]+)`\s*\|", line)
            if m:
                names.append(m.group(1))
        elif in_table:
            break
    return names


def config_fields():
    return {f.name for f in Config.__dataclass_fields__.values()}


# ----------------------------------------------------------------------------------------------------------
# Configuration options: Config dataclass <-> template.py <-> Simulation-configuration page
# ----------------------------------------------------------------------------------------------------------
def test_template_has_every_configuration_option():
    keys = set(re.findall(r"^(\w+)\s*=", read(ROOT / "configFiles" / "template.py"), flags=re.M))
    assert config_fields() - keys == set(), "options of Config missing in configFiles/template.py"
    assert keys - config_fields() == set(), "options of configFiles/template.py unknown to Config"


def test_configuration_page_documents_every_option():
    documented = set(first_column_names("Simulation-configuration.md", "The options of the configuration file are:"))
    assert config_fields() - documented == set(), "add these options to docs/wiki/Simulation-configuration.md"
    assert documented - config_fields() == set(), "options of docs/wiki/Simulation-configuration.md unknown to Config"


# ----------------------------------------------------------------------------------------------------------
# operadar() arguments <-> Execution page
# ----------------------------------------------------------------------------------------------------------
def test_execution_page_documents_every_argument_of_operadar():
    arguments = set(inspect.signature(operadar).parameters)
    documented = set(first_column_names("Execution-of-operadar.md", "# `operadar()` arguments"))
    assert arguments - documented == set(), "add these arguments to docs/wiki/Execution-of-operadar.md"
    assert documented - arguments == set(), "arguments of docs/wiki/Execution-of-operadar.md unknown to operadar()"


def test_examples_of_the_execution_page_use_existing_arguments():
    arguments = set(inspect.signature(operadar).parameters)
    code = "\n".join(re.findall(r"```python\n(.*?)```", read(WIKI / "Execution-of-operadar.md"), flags=re.S))
    for call in re.findall(r"\boperadar\(([^()]*)\)", code):
        used = set(re.findall(r"(\w+)\s*=(?!=)", call))
        assert used - arguments == set(), f"unknown arguments in an example of Execution-of-operadar.md: {call}"


def test_command_line_options_exist():
    cli = read(ROOT / "operadar" / "cli.py")
    for option in set(re.findall(r"(--\w+)", read(WIKI / "Execution-of-operadar.md"))) - {"--help"}:
        assert f'"{option}"' in cli or f"'{option}'" in cli, f"{option} is not an option of operadar/cli.py"


def test_example_configurations_are_listed():
    examples = {p.name for p in (ROOT / "configFiles").glob("test_*.py")}
    listed = set(re.findall(r"`(test_\w+\.py)`", read(WIKI / "Execution-of-operadar.md")))
    assert examples == listed, "update the table of example configurations in docs/wiki/Execution-of-operadar.md"


# ----------------------------------------------------------------------------------------------------------
# Files, scripts and links named in the documentation
# ----------------------------------------------------------------------------------------------------------
REPO_PATH = re.compile(r"(?<![\w/])(?:\./)?((?:configFiles|plot_tools|tables_generator|operadar|tests|docs)"
                       r"/[\w./-]*[\w/])(?![\w<-])")
# folders created by the user or at run time (not in the repository)
USER_FOLDERS = {"modelFiles", "operadarFiles", "operadar_env", "log", "libFortran"}
DOT_FOLDER = re.compile(r"(?<![\w.])\./([A-Za-z_][\w-]*)/")
SCRIPTS = re.compile(r"(?<![\w/.-])([\w-]+\.sh|operadar_multi\w*\.py)\b")


@pytest.mark.parametrize("page", DOCS, ids=lambda p: str(p.relative_to(ROOT)))
def test_files_named_in_the_documentation_exist(page):
    text = read(page)
    base = page.parent if page.parent != WIKI else ROOT
    missing = []
    for path in set(REPO_PATH.findall(text)) | set(SCRIPTS.findall(text)):
        if any(c in path for c in "<>*{}") or path.startswith("tables_generator/tables/") and path.count("/") > 3:
            continue                                         # placeholders and generated tables
        if not ((ROOT / path).exists() or (base / path).exists()
                or (path.startswith("tables_generator/tables/") and (ROOT / "tables_generator/tables").exists())):
            missing.append(path)
    for folder in set(DOT_FOLDER.findall(text)) - USER_FOLDERS:
        if not (ROOT / folder).is_dir() and not (base / folder).is_dir():
            missing.append(f"./{folder}/")
    for path in re.findall(re.escape(BLOB_URL) + r"([\w./-]+)", text):
        if not (ROOT / path).exists():
            missing.append(path)
    assert missing == [], f"{page.relative_to(ROOT)} names files that do not exist"


@pytest.mark.parametrize("page", DOCS, ids=lambda p: str(p.relative_to(ROOT)))
def test_links_to_wiki_pages_point_to_existing_pages(page):
    targets = re.findall(re.escape(WIKI_URL) + r"/([\w-]+)", read(page))
    missing = sorted({t for t in targets if not (WIKI / f"{t}.md").exists()})
    assert missing == [], f"{page.relative_to(ROOT)} links to wiki pages that are not in docs/wiki/"


def test_every_wiki_page_is_linked_from_another_page():
    text = "\n".join(read(p) for p in DOCS)
    orphans = [p.stem for p in PAGES if p.stem != "Home" and f"{WIKI_URL}/{p.stem}" not in text
               and p.stem not in ("FAQ", "ToDo-List")]
    assert orphans == [], "link these pages from another page (e.g. Home.md)"


# ----------------------------------------------------------------------------------------------------------
# Installation page <-> pyproject.toml
# ----------------------------------------------------------------------------------------------------------
def supported_pythons():
    spec = tomllib.loads(read(ROOT / "pyproject.toml"))["project"]["requires-python"]
    low = re.search(r">=\s*3\.(\d+)", spec)
    high = re.search(r"<\s*3\.(\d+)", spec)
    assert low and high, f"requires-python should have a lower and an upper bound: {spec}"
    return [f"3.{m}" for m in range(int(low.group(1)), int(high.group(1)))], f"3.{high.group(1)}"


def test_installation_page_gives_the_supported_python_versions():
    versions, first_unsupported = supported_pythons()
    text = read(WIKI / "Installation-tutorial.md")
    assert f"Python {' or '.join(versions)}" in text, f"Installation-tutorial.md should say Python {' or '.join(versions)}"
    assert f"Python {first_unsupported} is not supported" in text
