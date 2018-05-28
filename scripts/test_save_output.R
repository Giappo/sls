library(xlsx)
results_file      <- paste0(getwd(),"//results//table2.xlsx")
figure_name       <- paste0(getwd(),"//results//figura_errori.png")
figure.error_bars <- result$figure.error_bars
sheet_name        <- result$sheet_name

# create a new workbook for outputs
#++++++++++++++++++++++++++++++++++++
# possible values for type are : "xls" and "xlsx"
wb<-createWorkbook(type="xlsx")
# Define some cell styles
#++++++++++++++++++++++++++++++++++++
# Title and sub title styles
TITLE_STYLE <- CellStyle(wb)+ Font(wb,  heightInPoints=16,
                                   color="blue", isBold=TRUE, underline=1)
SUB_TITLE_STYLE <- CellStyle(wb) +
  Font(wb,  heightInPoints=14,
       isItalic=TRUE, isBold=FALSE)
# Styles for the data table row/column names
TABLE_ROWNAMES_STYLE <- CellStyle(wb) + Font(wb, isBold=TRUE)
TABLE_COLNAMES_STYLE <- CellStyle(wb) + Font(wb, isBold=TRUE) +
  Alignment(wrapText=TRUE, horizontal="ALIGN_CENTER") +
  Border(color="black", position=c("TOP", "BOTTOM"),
         pen=c("BORDER_THIN", "BORDER_THICK"))
# Create a new sheet in the workbook
#++++++++++++++++++++++++++++++++++++
sheet <- createSheet(wb, sheetName = sheet_name)
#++++++++++++++++++++++++
# Helper function to add titles
#++++++++++++++++++++++++
# - sheet : sheet object to contain the title
# - rowIndex : numeric value indicating the row to
#contain the title
# - title : the text to use as title
# - titleStyle : style object to use for title
xlsx.addTitle<-function(sheet, rowIndex, title, titleStyle){
  rows <-createRow(sheet,rowIndex=rowIndex)
  sheetTitle <-createCell(rows, colIndex=1)
  setCellValue(sheetTitle[[1,1]], title)
  setCellStyle(sheetTitle[[1,1]], titleStyle)
}
# Add title and sub title into a worksheet
#++++++++++++++++++++++++++++++++++++
# Add title
xlsx.addTitle(sheet, rowIndex=1, title="US State Facts",
              titleStyle = TITLE_STYLE)
# Add sub title
xlsx.addTitle(sheet, rowIndex=2,
              title="Data sets related to the 50 states of USA.",
              titleStyle = SUB_TITLE_STYLE)
# Add a table into a worksheet
#++++++++++++++++++++++++++++++++++++
d.f <- data.frame(result$results.table)
addDataFrame(d.f, sheet, startRow=3, startColumn=1,
             colnamesStyle = TABLE_COLNAMES_STYLE,
             rownamesStyle = TABLE_ROWNAMES_STYLE)
# Change column width
setColumnWidth(sheet, colIndex=c(1:ncol(state.x77)), colWidth=11)
# Add a plot into a worksheet
#++++++++++++++++++++++++++++++++++++
# create a png plot


png(figure_name, height=1000, width=1000, res=300, pointsize=8)
figure.error_bars
dev.off()
# Create a new sheet to contain the plot
# sheet <-createSheet(wb, sheetName = "boxplot")
# Add title
xlsx.addTitle(sheet, rowIndex=1, title="Box plot using InsectSprays data",
              titleStyle = TITLE_STYLE)
# Add the plot created previously
addPicture(figure_name, sheet, scale = 1, startRow = 4,
           startColumn = 1)
# remove the plot from the disk
res<-file.remove(figure_name)
# Save the workbook to a file...
#++++++++++++++++++++++++++++++++++++
saveWorkbook(wb, "table2")
