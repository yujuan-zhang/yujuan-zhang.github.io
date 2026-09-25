-- Computes real word count from the post body and inserts "X min read"
-- as the first block, right under the title/subtitle Quarto renders from
-- metadata. Applies to every post under blog/ via blog/_metadata.yml, so
-- it never has to be maintained by hand per-post.
local WORDS_PER_MINUTE = 200

function count_words(blocks)
  local words = 0
  local function count(el)
    if el.t == "Str" then
      words = words + 1
    end
  end
  for _, block in ipairs(blocks) do
    pandoc.walk_block(block, { Str = count })
  end
  return words
end

function Pandoc(doc)
  local words = count_words(doc.blocks)
  local minutes = math.max(1, math.ceil(words / WORDS_PER_MINUTE))
  local note = pandoc.Div(
    { pandoc.Plain({ pandoc.Str(minutes .. " min read") }) },
    pandoc.Attr("", { "reading-time" })
  )
  table.insert(doc.blocks, 1, note)
  return doc
end
