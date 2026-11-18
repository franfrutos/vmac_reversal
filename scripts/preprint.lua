-- Figures of the preprint (journal mode): all of them span both columns, centered, with
-- these widths (of the page) so that the tallest ones do not take a whole page

local widths = {
  ["fig-1"] = "66%",
  ["fig-2"] = "92%",
  ["fig-3"] = "68%",
  ["fig-4"] = "66%",
  ["fig-5"] = "82%",
  ["fig-6"] = "66%",
  ["fig-7"] = "82%"
}

function Div(div)
  if not div.identifier:match("^fig%-") then
    return nil
  end

  div.attributes["apa-twocolumn"] = "true"
  return div:walk {
    Image = function(img)
      img.attributes["apa-twocolumn"] = "true"
      img.attributes["fig-pos"] = nil -- [H] does not work with figures across columns
      img.attributes["width"] = widths[div.identifier] or img.attributes["width"]
      img.attributes["fig-align"] = "center"
      return img
    end
  }
end
