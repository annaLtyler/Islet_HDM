#adds a commented out header to a data table.

write_table_with_header <- function(data.table, file.name, header = "",
  sep = "\t", col.names = TRUE, row.names = FALSE){
  write.table(data.table, file.name, quote = FALSE, sep = sep, 
    col.names = col.names, row.names = row.names)
  fConn <- file(file.name, 'r+')
  Lines <- readLines(fConn)
  writeLines(c(header, Lines), fConn)
  close(fConn)
}
