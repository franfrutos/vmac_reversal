-- -----------------------------------------------------------------------------
-- Author: Francisco Garre-Frutos
-- -----------------------------------------------------------------------------

-- Word count (text, notes and abstract) on the title page, with the
-- quarto-wordcount extension (Andrew Heiss). Its filters run on a copy of the
-- document, because apaquarto runs citeproc again and the references would be
-- duplicated

local ext_dir = pandoc.path.join({
  pandoc.path.directory(PANDOC_SCRIPT_FILE), "..", "_extensions", "andrewheiss", "wordcount"
})

-- Load one of the extension's filters in its own environment and return its Pandoc function
local function load_filter(file)
  local env = setmetatable({}, { __index = _G })
  local chunk = assert(loadfile(pandoc.path.join({ ext_dir, file }), "t", env))
  chunk()
  return env.Pandoc
end

local function num(x)
  if x == nil then return 0 end
  if type(x) == "number" then return x end
  return tonumber(pandoc.utils.stringify(x)) or 0
end

local function with_commas(n)
  local s = tostring(math.floor(n))
  local out = s:reverse():gsub("(%d%d%d)", "%1,"):reverse()
  return (out:gsub("^,", ""))
end

-- Figures and tables are still Divs (#fig-..., #tbl-...) at this stage, so drop them
-- before counting to leave their captions out
local function drop_floats(blocks)
  return pandoc.walk_block(pandoc.Div(blocks), {
    Div = function(d)
      if d.identifier:match("^fig%-") or d.identifier:match("^tbl%-") then
        return {}
      end
    end
  }).content
end

function Pandoc(doc)
  local ok, err = pcall(function()
    local run_citeproc = load_filter("citeproc.lua")
    local run_count = load_filter("wordcount.lua")
    local copy = pandoc.Pandoc(drop_floats(doc.blocks:clone()), doc.meta)
    copy = run_count(run_citeproc(copy))
    local n = num(copy.meta.wordcount_body_words)
      + num(copy.meta.wordcount_note_words)
      + num(copy.meta.wordcount_abstract_words)
    if n > 0 then
      doc.meta.note = pandoc.Inlines({ pandoc.Str("Word count: " .. with_commas(n)) })
    end
  end)
  if not ok then
    io.stderr:write("wordcount_titlepage.lua: word count not added (" .. tostring(err) .. ")\n")
  end
  return doc
end
