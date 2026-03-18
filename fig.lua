-- shortcap-subfig.lua
-- Allows use of fig-scap with figures containing subfigures (fig-subcap)

local is_latex = FORMAT:match("latex")

function Div(el)
  if not is_latex then
    return nil
  end
  if not el.identifier or el.identifier == "" then
    return nil
  end

  -- Only act on figures with ids starting with "fig-"
  if not el.identifier:match("^fig%-") then
    return nil
  end

  -- The final paragraph inside the Div is interpreted as the long caption
  local blocks = el.content
  local longcap = nil
  if #blocks > 0 then
    local last = blocks[#blocks]
    if last.t == "Para" or last.t == "Plain" then
      longcap = pandoc.utils.stringify(last)
      table.remove(blocks)
    end
  end

  local shortcap = el.attributes["fig-scap"]

  local lines = {}
  table.insert(lines, "\\begin{figure}")
  table.insert(lines, "\\centering")

  -- Render the regular content (subfigures) in LaTeX
  for _, block in ipairs(blocks) do
    table.insert(lines, pandoc.write(pandoc.Pandoc({ block }), "latex"))
  end

  if longcap then
    if shortcap then
      table.insert(lines, string.format("\\caption[%s]{%s}", shortcap, longcap))
    else
      table.insert(lines, string.format("\\caption{%s}", longcap))
    end
  end

  table.insert(lines, string.format("\\label{%s}", el.identifier))
  table.insert(lines, "\\end{figure}")

  return pandoc.RawBlock("latex", table.concat(lines, "\n"))
end
