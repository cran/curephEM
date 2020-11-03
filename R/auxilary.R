check.binary = function(col)
{
  tmp = unique(col)
  return(sum(!is.na(tmp))<=2)
}
