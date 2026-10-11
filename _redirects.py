# Post-render: blog folders were renamed from 260924my-story to
# 2026-09-24my-story. Write small redirect pages at the old paths so
# previously shared links keep working.
import os, re

SITE = "_site/blog"
STUB = ('<!doctype html><html><head><meta charset="utf-8">'
        '<title>Redirecting</title><link rel="canonical" href="{url}">'
        '<meta http-equiv="refresh" content="0; url={url}"></head>'
        '<body><a href="{url}">{url}</a></body></html>\n')

for new in os.listdir(SITE):
    m = re.fullmatch(r"20(\d\d)-(\d\d)-(\d\d)(.+)", new)
    if not m or m.group(1) + m.group(2) + m.group(3) > "261009":
        continue
    old = os.path.join(SITE, "".join(m.groups()))
    os.makedirs(old, exist_ok=True)
    for page in ("index.html", "index-chinese.html"):
        if os.path.exists(os.path.join(SITE, new, page)):
            with open(os.path.join(old, page), "w") as f:
                f.write(STUB.format(url="/blog/%s/%s" % (new, page)))
