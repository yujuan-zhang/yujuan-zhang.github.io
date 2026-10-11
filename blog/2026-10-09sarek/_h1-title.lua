-- These posts have no front matter: use the leading "# ..." heading as the
-- page title (rendered in the site's title block) instead of editing the file.
-- Quarto fills in the file name as a fallback title, and other filters may
-- insert blocks ahead of the heading, so look for the first header instead.
function Pandoc(doc)
  for i, b in ipairs(doc.blocks) do
    if b.t == "Header" then
      if b.level == 1 then
        doc.meta.title = pandoc.MetaInlines(b.content)
        doc.meta.pagetitle = pandoc.utils.stringify(b.content)
        doc.blocks:remove(i)
      end
      break
    end
  end
  return doc
end
