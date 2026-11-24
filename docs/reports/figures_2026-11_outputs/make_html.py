"""Self-contained HTML version (figures embedded) of a Markdown report of docs/reports/.

    python docs/reports/figures_2026-11_outputs/make_html.py docs/reports/2026-11_operadar_output_analysis.md

Needs the `markdown` package (pip install markdown). The HTML files are not tracked by git.
"""
import base64
import os
import re
import sys

import markdown

CSS = """body{font-family:system-ui,sans-serif;max-width:1150px;margin:2em auto;padding:0 1em;line-height:1.5;color:#222}
img{max-width:100%;border:1px solid #ddd;margin:.5em 0}
table{border-collapse:collapse;margin:1em 0;font-size:.92em}
th,td{border:1px solid #ccc;padding:.3em .6em;vertical-align:top}
th{background:#f3f3f3}
code{background:#f4f4f4;padding:0 .2em;border-radius:3px}
pre{background:#f4f4f4;padding:.8em;overflow-x:auto}
h1{border-bottom:2px solid #2a6fdb}
h2{border-bottom:1px solid #ccc;margin-top:2em}"""


def main(md_path):
    folder = os.path.dirname(os.path.abspath(md_path))
    text = open(md_path, encoding='utf-8').read()
    title = text.splitlines()[0].lstrip('# ').strip()
    body = markdown.markdown(text, extensions=['tables', 'fenced_code', 'toc'])

    def embed(m):
        src = m.group(1)
        path = os.path.join(folder, src)
        if not os.path.isfile(path):
            print(f'/!\\ missing figure {src}')
            return m.group(0)
        data = base64.b64encode(open(path, 'rb').read()).decode()
        return f'src="data:image/png;base64,{data}"'

    body = re.sub(r'src="([^"]+\.png)"', embed, body)
    html = (f'<!DOCTYPE html><html lang="en"><head><meta charset="utf-8"><title>{title}</title>'
            f'<style>\n{CSS}\n</style></head><body>{body}</body></html>\n')
    out = os.path.splitext(md_path)[0] + '.html'
    open(out, 'w', encoding='utf-8').write(html)
    print(f'{out}: {os.path.getsize(out) / 1e6:.1f} MB')


if __name__ == '__main__':
    main(sys.argv[1])
