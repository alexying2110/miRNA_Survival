###Libraries###
#Used for llply, faster version of lapply
library(plyr)
#Used for full_join, which is slightly slower than the match join
library(dplyr)
#Used for the JIT compiler
library(compiler)
#Used for data tables, keys, fread, and fwrite
library(data.table)
#Used for all survival analysis functions
library(survival)
#Used for nicer looking plots
library(ggplot2)
#Used to allow ggplot to work in conjunction with survival
library(ggfortify)


###Set Up###
#Turns on JIT compiler which makes the whole thing a bit faster, this is the most compiled
#setting for the JIT
enableJIT(3)


###Helper Functions###
#Makes the data table from the given folder name using fread, which is significantly faster
#than read.table. Subsets the data by only pulling the columns needed (pulled by name to
#increase security). Renames the column based on the file name. Returns only the expression
#column, but first aligns the data based where the rows match up with the vector "miRNA"
#(using the match function).An empty bracket is at the end of the return line because
#the data.table package is slightly bugged and this allows for proper returning
#Unfortunately, while this method is significantly faster than full_join, it is closer to
#a left join, so miRNAs that appear in vector of miRNA but don't appear in the file will
#be set to NA, but miRNAs that appear in the file, but not in the vector will remain
#unassigned.
#miRNA is a vector of miRNA names
#file is a file name
getExpression <- function(miRNA, file){
	dt <- fread(file, sep = "\t")[,c("miRNA_ID", "reads_per_million_miRNA_mapped")]
	names(dt)[2] <- file
	dt[na.omit(match(miRNA,dt$miRNA_ID)), 2][]
}

#makes a data table like the function above except returns a two column data frame with the
#second column renamed to be the file name
#file is a file name
makeDT <- function(file){
	dt <- fread(file, sep = "\t")[,c("miRNA_ID", "reads_per_million_miRNA_mapped")]
	names(dt)[2] <- file
	dt
}

#Produces a line of values for the statistical analysis of the data. Creates a survdiff
#object, which is processed to determine how distinct the high expression and low
#expression survival curves are. The function also produces a cox object, summarizes it
#and extracts the relevant info.
#DT is a data table (this argument isn't used in normal operation)
#surv is the surv object created in the parten function
#line is the binarized vector of the same length as surv
survivalSummary <- function(DT, surv, line){
	#subset and na.action arguments are pretty much useless, it just gets rid of a warning
	#in RStudio that annoys me. Cofactors can be added into the formula but aren't necessarily
	#productive with such a low sample size
	sdObj <- survdiff(formula = surv~line,subset = T, na.action = "na.omit")

	#Return vector, returns cox Hazard Ratio, cox p value (probability that the hazard ratio is 1),
	#and the survival difference p value (probability that the curves are not distinct)
	c(summary(coxph(surv~line))$coefficients[c(2,5)], 1 - pchisq(sdObj$chisq, length(sdObj$n) - 1))
}

#Produces a line of values for the statistical difference between the two sections of the
#data (high expressions and low expressions). This function also returns the median expression
#the mean of the lower and higher expressions, and the reciprocal of the fold, as well as the
#adjusted p value
#exp is the expression data
#binary is the previously calculated binarized data
#Both values are a vector of the same length and are intended to be used in mapply
tTests <- function(exp, binary){
	#Splits exp into a list with position determined by the binarized data
	splitExp <- split(exp, binary)
	#Two sample t-test using strata split above. This is significantly faster than using formula
	#notation in t.test. 0 and 1 are named specifically to ensure the order is correct for the
	#t-test.
	tResult <- t.test(splitExp$'0', splitExp$'1')
	#Extracts the means from tResult and casts to numeric
	means <- as.numeric(tResult$estimate)
	#Extracts the p-value from tResult
	pVal <- tResult$p.value
	#Returns the median expression, low mean, high mean, reciprocal of the fold (to prevent
	#divide-by-0 errors, the p-value, and the fdr adjusted p-value
	c(median(exp), means[1], means[2], means[1]/means[2], pVal, p.adjust(pVal, "fdr") )
}

###TODO: Finish up this method to clean up the plot and make it more presentable with the title
###as the miRNA name
savePDFs <- function(survObj, binary){
	SF <- survfit(survObj ~ binary)
	names(SF$strata) <- c("Low Expression", "High Expression")
	autoplot(SF, conf.int = F, surv.color = c(132, "red"))
}


###TODO: PUT THESE FUNCTIONS DIRECTLY INTO THE CODE
#Although this function uses a for loop, it is memory and time efficient to a reasonable degree
#as it modifies the data table in place, removing the primary weakness of a for loop in R.
#The if statement verifies that the modifications aren't attempted on a data table without
#duplicates
#DT is the data.table with duplicate columns
#Note: this is maybe not the single fastest way to accomplish the goal, but it executes
#pretty much instantly and there's no point to optimize further at the cost of more memory
dupCollapse <- function(DT){
	dups <- unique(names(DT)[duplicated(names(DT))])
	if (length(dups) > 0) {
		for (i in 1:length(dups)) {
			dupVec <- names(DT) == dups[i]
			DT[, which(dupVec)[1] := rowMeans(DT[, which(dupVec), with = F])]
			DT[, which(dupVec)[-1]  := NULL]
		}
	}
	DT[]
}

#Transforms the input data table into a data table that has expression values replaced with binary values representing
#whether the expression value is greater than or less than the median expression level, respectively. The patient id's
#are reintroduced in the same order as the data table was split, and each column is assigned an miRNA_id.
#llplying the split matrix is the fastest way to operate on the each row, followed by do.calling cbind over the list.
#Converting to data table has very low overhead, and the two operations on the data table are passed by reference.
#Memory usage is under .2 mB.
#DT is the cleaned data table
binarize <- function(DT){
	binaryDT <- as.data.table(do.call(what = "cbind",args = ,llply(split(as.matrix(DT[,-1]), f = row(DT[,-1])), function(x) if_else(x > median(x), 1, 0))))
	setnames(binaryDT, names(binaryDT), unlist(DT[,1]))
	binaryDT[, "patient_id" := names(DT)[-1]][]
}


###Main Functions###
#Execute this to create the data table with all the patient data. Fast makes the data table 4 times faster without it,
#result is the same regardless but setting fast to F uses a method that is known to work without mistakes and uses no
#short cuts, and still operates about 20 times faster. Weirdly enough when using all.equal() on the result of this
#function to test against the old function, the result is FALSE, but checking the individual values is correct, in any
#case, the results remain the same so it is likely just a type bug somewhere in the old code. Setting fast = F results
#in 34% the memory usage (inconsistent based on how many times the call has been made previously)
#path is the folder where the subfolders are located.
#meta_data is the
makeExpTable <- function(path = ".", meta_data, fast = T, csv = F){
	#Sets wd to where all the files are. Code only functions if the the folder selected contains
	#the subfolders that each contain the text file with the miRNA data
	setwd(path)
	#The metadata of the files, used for the patient_id information
	#Assumed to be a csv, if not, remove the sep statement
	metaData <- fread(meta_data, sep = ",")
	#List all files that match pattern in the wd. The patter is selected so that only text files
	#with miRNA data are pulled.
	fileList <- unlist(list.files(path = ".", pattern = "*mirnas.quantification.txt$", recursive = T))
	if (fast) {
		#Creates a base data frame that all other data sets will be compared to, no subsetting performed as
		#it will unnecessarily increase operation time by a few miliseconds
		miRNAs <- fread(fileList[[1]])[,"miRNA_ID"]
		#First it llply's the file list into a list of expression data values that line up with the miRNA
		#list extracted above. llply() is significantly faster than lapply. This list of expression data values
		#in the form of data table columns is not saved into its own variable. It's assumed that this will
		#require less memory. Tested on the LUSC data set, even if the list was saved, it would be a trivial
		#amount of megatbytes. The miRNA list and data table list are joined to forma  single list, which is converted
		#to a data frame. Since the alignment happens in the creation of the list, no additional joining is needed,
		#which greatly increases speed. getExpression() increases operation speed by reducing miRNAs to a vector
		expDT <- as.data.table(c(miRNAs,llply(fileList, function(x) getExpression(getElement(miRNAs, "miRNA_ID"), x))))
	} else {
		#Very similar to above except uses Reduce() to recursively apply the full_join on the list of data frames
		#created by the llply call. Returns a final data table with all the columns merged together.
		expDT <- setDT(Reduce(function(x, y) full_join(x, y, by = "miRNA_ID"), llply(fileList, makeDT)))
	}
	#Renames the columns based on the patient id's in the metaData.
	#First takes all the file names, splits them by the "/", uses vapply (faster version of sapply, but have to
	#specify return format) to grab just the folder name, finds where the folder name matches up with the
	#corresponding column in the meta data, uses that to get a row number, then uses matrix subsetting to extract
	#the column in the meta data in order of the folder names. Vectorizes then sets it all to lowercase.
	###AS OF AN UPDATE TO THE TCGA METADATA FILES IT APPEARS THAT THE COLUMNS HAVE BEEN RENAMED AND THIS LINE NO LONGER WORKS
	#setnames(expDT, names(expDT)[2:ncol(expDT)] , tolower(metaData$"cases/0/samples/0/submitter_id"[match(vapply(strsplit(names(expDT)[2:ncol(expDT)], "/"), "[", 1, FUN.VALUE = ""), metaData$file_id)]))
	
	#Updated renamer 
	setnames(expDT, names(expDT)[2:ncol(expDT)] , tolower(metaData$"associated_entities/0/entity_submitter_id"[match(vapply(strsplit(names(expDT)[2:ncol(expDT)], "/"), "[", 1, FUN.VALUE = ""), metaData$file_id)]))
	
	#If the user choses to save, fwrite writes the csv to the same folder as the data folder and names it.
	if (csv) {
		fwrite(expDT, "miRNA_Expression_Data.csv")
	}
	expDT
}

#First converts the datax table to a data frame, which makes extracting columns slightly faster, this
#conversion isdone by reference to save memory. The colnames are subset and "-01" is searched for(subsetting
#makes the process marginally faster), and those columns are extracted. miRNA that have an average expression
#less than 1 get removed. The data frame is converted back to data table and duplicate columns are averaged.
#Submitter id is removed.
#DT is the data table produced by makeExpTable()
rmLowExp <- function(DT){
	#Converting and unconverting saves 60 ms on a 1881 by 568 data table.
	setDF(DT)
	#Does not use the data table := operator to remove columns because it's significantly slower.
	#Extract only columns with "-01"
	DT <- DT[which(rowMeans(DT[,-1]) > 1),c(1,grep("-01",substr(names(DT), 13, 16)))]
	#empty bracket to return properly
	setDT(DT)
	#Just removes extraneous tag that is no longer needed for survival operations
	setnames(DT, names(DT)[-1], substr(names(DT)[-1], 9, 12))
	dupCollapse(DT)
	DT[]
}

#Takes data table of expression values to produce a table of survival statistics and statistics describing the difference
#between the high expression and low expression levels.
#DT is the cleaned data.table
#clinicalCSV is the file location of the clinical data
#savePDF determines whether the significant survival curve plots are saved to the same directory the csv is located in a
#folder titled "Survival_Plots"
survivalProcess <- function(DT, clinicalCSV, savePDF = F) {
	#Quickly read given file, assume separators but may need a modification if not actually the case
	clinDT <- fread(clinicalCSV, sep = ",")
	#Reduce clinDT to only the patients with expression data in DT
	clinDT <- clinDT[patient_id %in% names(DT), ]
	#Changes the column order of DT to match the order of patient id's in clinDT, this allows each expression data set to
	#be processed quickly outside of the data column str ucture while retaining the ability to use clinDT since the orders
	#are identical
	setcolorder(DT, c(1, match(clinDT$patient_id, names(DT))))

	#Produces the Surv object for recurrence events and overall survival events
	survObsRecurr <- Surv(clinDT$survival_recurrence, clinDT$recurrence)
	survObsDeath <- Surv(clinDT$survival_overall, clinDT$live_death)

	#Subsets the data.table to ignore the first column, which is a character vector, to allow the entire data table to quickly
	#be cast to matrix (this makes splitting much faster). Splits the table by rows to produce a vector of expression values for
	#each miRNA and produces a list
	expList <- split(as.matrix(DT[,-1]), f = row(DT[,-1]))
	#Goes through the expression value list and binarizes using the fast llply function
	binaryList <- llply(expList, function(x) if_else(x > median(x), 1, 0))

	#First uses mapply to apply a combination function of survivalSummary and tTests across the binarized list and the
	#expression data list. Mapply is used because it is significantly faster to use t.test with two samples than to use a formula
	#(see comments in tTests function). c is used to glue all the function outputs together, and the result is rounded. mapply
	#produces a matrix with each column representings the function output, which is immediately transposed and cast to data table.
	statsDT <- as.data.table(t(mapply(function(x,y) round(c(survivalSummary(clinDT, survObsRecurr, y), survivalSummary(clinDT, survObsDeath, y), tTests(x, y)),4), expList, binaryList)))

	#Assigns column names, see functions for what is what
	setnames(statsDT, names(statsDT), c("Recurr.HR", "Recurr.p", "Recurr.SD.p", "Death.HR", "Death.p", "Death.SD.p", "med.exp", "low.mean", "high.mean", "1/fold", "p.val", "fdr"))
	#Creates a column for the miRNA_ID's
	statsDT[, "miRNA_ID" := DT$miRNA_ID]
	#Moves the miRNA_ID column to the front
	setcolorder(statsDT, c(ncol(statsDT), 1:(ncol(statsDT) - 1)))
		
	#Saves the pdfs of significant plots
	if (savePDF) {
		sigIndex <- which(statsDT$Recurr.SD.p > .05)
		goodMiRNA <- statsDT[sigIndex, "miRNA_ID"]
		l_ply()
	}


	#Typical return line
	statsDT[]
}

survivalProcess2 <- function(DT, clinicalCSV){


}


#Returns a list of miRNA worth looking at the literature for, listed with what value is significant.
miRNAForInvestigation = function(DT){
	list("Recurrence Survival Difference" = DT$miRNA_ID[which(DT$Recurr.SD.p <= .05)], 
		 "Recurrence Cox" = DT$miRNA_ID[which(DT$Recurr.p <= .05)],
		 "Death Survival Difference" = DT$miRNA_ID[which(DT$Death.SD.p <= .05)], 
		 "Death Cox" = DT$miRNA_ID[which(DT$Death.p <= .05)])
}




