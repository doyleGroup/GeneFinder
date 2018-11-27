## Welcome to GeneFinder

GeneFinder is a python program enabling fast discovery of proteins, noncoding RNAs, and phenotypes from user-provided article databases. Gene Finder was written by Samuel Berkun, Andrew Childers, and Adele Doyle.

The complete text of GeneFinder.py is shown below. 

```markdown
`# Language: Python 3
# GeneFinder.py
# Authored by Samuel Berkun and Andrew Childers of the Doyle research group at the University of California Santa Barbara
# Purpose: To compare an excel sheet of keywords (gene symbols and names) with another
# excel sheet of text (article titles and abstracts) in order to return
# matches of keywords found and where they were found. Genes and scientific
# papers do not need to be the input data, however the code is optimized
# for the nuances associated with searching for genes. The resulting gene
# counts are quantified and genes are analyzed by running the GeneAnalysis 
# function at or after the conclusion of this script. This script yields a
# list of the indices of genes found and articles they are found in. 
# Copyright: GNU AGPLv3
# Permissions of this strongest copyleft license are conditioned on making available complete source code of licensed works and modifications, which include larger works using a licensed work, under the same license. Copyright and license notices must be preserved. Contributors provide an express grant of patent rights. When a modified version is used to provide a service over a network, the complete source code of the modified version must be made available.



import pandas as pd   #excel file reading
import xlrd           #in case pandas didn't import it for some reason
import xlsxwriter     #in case pandas didn't import it for some reason
import numpy as np    #efficient boolean matrix
import time           #internal timer
import re             #regex
import tkinter        
import tkinter.ttk    #GUI
import matplotlib.pyplot as plt; 
from matplotlib.ticker import FuncFormatter



def GeneFinder():
    
    #makes a progress bar that updates while the function is running
    #tkinter needs to be imported
    def trackProgress(bigFunction,totalLength):
        #bigFunction is a function that takes 1 parameter, updateProgressBar
        #totalLength is an int
        try:
            ROOT = tkinter.Tk();
            Lwidget = tkinter.Label(ROOT,text="");
            Pwidget = tkinter.ttk.Progressbar(ROOT,orient="horizontal",
                                         length=500, mode="determinate");
            Lwidget.pack();
            Pwidget.pack();
            Pwidget["maximum"] = totalLength;
            def updateProgressBar(num,txt):
                #every so often, the function calls updateProgressBar
                #to show the user how much progress it made
                Lwidget["text"] = txt; 
                Pwidget["value"] = num;
                ROOT.update();
                if num==totalLength: 
                    ROOT.destroy();
            def start():
                bigFunction(updateProgressBar);
            updateProgressBar(0,"");
            ROOT.after(1,start); #wait 1 ms so the window opens before start
            ROOT.mainloop();
        except:
            #if graphics don't work, we just print every once in a while
            def printProgress(num,txt): 
                if num==int(totalLength*0.25): print("25% done:  "+txt);
                if num==int(totalLength*0.50): print("50% done:  "+txt);
                if num==int(totalLength*0.75): print("75% done:  "+txt);
                if num==totalLength:           print("100% done: "+txt);
            bigFunction(printProgress);
            
    
    
    
    #A bad word is a gene symbol that also happens to be a common word.
    #for example, AS represents the gene ankylosing spondalitus
    #because they are common words, the program makes does not detect them
    def isBadWord(word): 
        badwords = [
            "AS","II","KIT","MICE","WARS","DO","ON","IMPACT",
            "WAS","CO","IV","P","MS","TH","MED","PER","ED","NM",
            "ATM","SET","MET","CAPS","SD","CA1","CA2","MSC","GOV"
        ]
        return word in badwords;
    
    
    #binary search
    #returns the index of the gene if found, or -1 if not
    def search(geneList,name):
        a = 0; b = len(geneList)-1;
        r = int((a+b)/2); #floor division
        while a<=b:
            if geneList[r]==name: return r;
            if geneList[r]>name:
                b = r-1;
            else:
                a = r+1;
            r = int((a+b)/2); #floor division
        return -1;
    
    
    #returns the index of the correct gene
    #see documentation
    def findGeneidx(wordTable, firstWord, geneList, stTitle):
        #in this section, we want the range of all the genes that could start
        #with firstWord
        rm = search(wordTable[0], firstWord);
        x = rm; #Lower bound
        y = rm; #Upper bound
        if rm<0: return rm;
        while(x>0 and wordTable[0][x-1]==firstWord): 
            x = x-1;
        while(y+1<len(wordTable[0]) and wordTable[0][y+1]==firstWord): 
            y = y+1;
        
        #we now use the range of genes that start with firstWord, and select the
        #longest match
        for c in range(y,x-1,-1):  #backwards loop allows us to select the last (longest) one
            geneidx = wordTable[1][c];
            if stTitle.startswith(geneList[geneidx]+' '): 
                return geneidx;
        
        #if we are here, that means that there exists a gene/genes that start
        #with firstWord, but none of them are matches
        return -1; 
    
    
    #does a gene share its symbol with another gene?
    #true if it shares a gene symbol, false otherwise
    def needsFlag(geneTable, idx):
        return ((idx>0 and geneTable[0][idx]==geneTable[0][idx-1]) or 
            (idx+1<len(geneTable[0]) and geneTable[0][idx]==geneTable[0][idx-1]))
    
    ###############################################################################
        
    
    
    print('Program starting!')
    #Variables generated by this code 
    #(aside from ArticleData, GeneData, matchidxList, and output) 
    #are automatically cleared at the end of the program as to not crowd the 
    #workspace and use up memory needed for the code to run quickly.
    
    time1 = time.time(); #setting up interal timer 
    
    
    print('This script assumes you have two .xlsx files:')
    print('1 Excel file with articles titles or titles+abstracts and 1 with the gene names.')
    print('Please add .xlsx to the end of every file name you type in.'
          +' Do not put quotes around any of your input.')
    convert = input('Would you like to process abstracts?'+
                    ' Type 1 for titles+abstracts, 0 for just titles. ')=='1';
    writeAfter = input('Would you like to save results to a spreadsheet?'+
                    ' Type 1 for yes, 0 for no. ')=='1';
    if(writeAfter):
        resultsfilename = input('What would you like the file to be named? '+
                                'If you type the name of an existing file, '+
                                'the results will be saved to that file. '+
                                'If you type the name of an existing file, '+
                                'the file must not be open on your computer. ');
                                
    filterGeneCat = input('Would you like to filter genes by locus group?'+
                    ' Type 1 for yes, 0 for no. ')=='1';
    if filterGeneCat:
        locusSelected = [False,False,False,False,False,False];
        print('protein-coding gene = 1');
        print('non-coding RNA = 2');
        print('phenotype = 3');
        print('pseudogene = 4');
        print('withdrawn = 5')
        print('other = 6');
        locuses = input('Type the numbers of the groups you need (Example: 126): ');
        for a in range(6):
            locusSelected[a] = locuses.count(str(a+1))>0;
    # Make sure to have the excel files in the template format prior to
    # entering in the file names next.
    articlefilename=input('Enter name of file containing article titles. ');
        # C:/Users/Cheese/Documents/MATLAB/A2.xlsx
    genefilename=input('Enter name of file containing gene names. '); 
        # C:/Users/Cheese/Downloads/BiologicalDataSet.xlsx
    # Load spreadsheets
    print('...');
    time2 = time.time();
    print('Time elapsed during user input: '+str(time2-time1)+' seconds.');
    
    
    ###############################################################################
    
    
    
    print('Fetching files...');
    dfG = pd.ExcelFile(genefilename); dfG = dfG.parse(dfG.sheet_names[0]);
    dfA = pd.ExcelFile(articlefilename); dfA = dfA.parse(dfA.sheet_names[0]);
    #dfA and dfG are now pandas dataframes, we want them as lists
    
    GeneData = [dfG[dfG.columns[0]].map(str),
                dfG[dfG.columns[1]].map(str),
                dfG[dfG.columns[2]].map(str),
                dfG[dfG.columns[3]].map(str)];
    ArticleData = ([dfA[dfA.columns[0]].map(str),
                    dfA[dfA.columns[1]].map(str)] 
                if convert else dfA[dfA.columns[0]].map(str));
    
    
    print('Preprocessing input...');
    
    def cleanPunctuation(ar):
        #this function does 3 things:
        # 1. removes punctuation
        # 2. removes extra spaces, as the punctuation was replaced with spaces. 
        #    This is important, because many gene names are multiple words
        # 3. trims the string. This is for genes, so they don't get false negatives
        #    from trying to match whitespace 
        return [str.strip(re.sub('  +',' ',re.sub('[^a-zA-Z0-9]',' ',e))) for e in ar];
    
    if convert:
        #clean the columns, then contenate them
        #concatenate is after so we can use a double space to keep titles and abstracts seperate
        ArticleData[0] = cleanPunctuation(ArticleData[0]);
        ArticleData[1] = cleanPunctuation(ArticleData[1]);
        ArticleData = [a+"  "+b for a,b in zip(ArticleData[0],ArticleData[1])];
    else:
        ArticleData = cleanPunctuation(ArticleData);
        
    ArticleData = list(set(ArticleData)); # unique(ArticleData)
    ArticleData = sorted(ArticleData);    # so we can easily find correct article in results
    
    def sortByFirstRow(arr):
        #backwards loop so we sort first row last, since first row sets order
        for a in range(len(arr)-1,-1,-1): 
            arr[a] = [x for _,x in sorted(zip(arr[0],arr[a]))];
        return arr;
    
    
    GeneData[0] = [elem.replace('~withdrawn','') for elem in GeneData[0]];
    for a in range(4):
        GeneData[a] = GeneData[a][1:]; #first row is just more column headers
        GeneData[a] = cleanPunctuation(GeneData[a]);
        GeneData[a] = list(map(str.upper, GeneData[a]));
        #the articles get uppercased too, but inside the loop
    
    if(filterGeneCat): #if they wanted to filter by locus
        keepGenes = [(locusSelected[0] and a=="PROTEIN CODING GENE" or
                      locusSelected[1] and a=="NON CODING RNA" or
                      locusSelected[2] and a=="PHENOTYPE" or
                      locusSelected[3] and a=="PSEUDOGENE" or
                      locusSelected[4] and a=="WITHDRAWN" or
                      locusSelected[5] and a=="OTHER") for a in GeneData[3]];
        for a in range(4):
            GeneData[a] = [elem for b,elem in enumerate(GeneData[a]) if keepGenes[b]]; 
    
    GeneData = sortByFirstRow(GeneData);
    
    #geneSymbFW and geneNameFW store the first word of every gene symbol or gene name.
    #this allows for efficient checking 
    #the second column is the index in GeneData that the gene corresponds to
    geneSymbFW = [GeneData[0],list(range(len(GeneData[0])))];
            # gene symbols, and the index of each gene symbol
    geneNameFW = [[a for a,b in zip(GeneData[1],GeneData[2]) if b=="APPROVED"],
                [a for a,b in       enumerate(GeneData[2]) if b=="APPROVED"]];
                 
            #only keep gene names from approved genes, since for withdrawn genes
            #the spreadsheet says "entry withdrawn" instead of the gene name
    
    geneSymbFW = sortByFirstRow(geneSymbFW); # so binary sort can be performed
    geneNameFW = sortByFirstRow(geneNameFW);
    
    geneSymbFW[0] = [a.partition(' ')[0] for a in geneSymbFW[0]]; #first word of each
    geneNameFW[0] = [a.partition(' ')[0] for a in geneNameFW[0]]; #first word of each
    
    
    time3 = time.time();
    print("Time elapsed during fetching and preprocessing: "+str(time3-time2)+" seconds.");
    
    
    
    ###############################################################################
    
    
    
    print("Finding articles containing 1 or more gene names...");
    
    
    matchidxList = [];
    
    def findGenes(updateProgressBar):
        empties = np.full(len(GeneData[0]),False);
        p=0;
        for a in range(len(ArticleData)): #for every article:
            #p will help us keep track of how many results we get for this article
            p=len(matchidxList); 
            
            hasSymb = empties.copy(); #whether the article has a given gene symbol
            hasName = empties.copy(); #whether the article has a given gene name
            
            articleText = ArticleData[a].upper()+" ";
            while(len(articleText)>0):
                #loop through each word of the article
                firstWord = articleText[:articleText.find(" ")];
                if firstWord=="": 
                    articleText = articleText[1:];
                    continue;
                    
                    
                #see if there is a gene name at this word
                geneidx = findGeneidx(geneNameFW, firstWord, GeneData[1], articleText);
                if geneidx>0: 
                    #we found a gene name at this word
                    #if it's a new gene for this article, add it to the list of matches
                    if not (hasSymb[geneidx] or hasName[geneidx]): 
                        matchidxList.append(geneidx);
                    hasName[geneidx] = True;
                    #go to next word in the article
                    articleText = articleText[len(GeneData[1][geneidx])+1:];
                    continue;
                    
                #see if there is a gene symbol at this word
                if isBadWord(firstWord): #badwords would be false positives
                    articleText = articleText[len(firstWord)+1:];
                    continue;
                geneidx = findGeneidx(geneSymbFW, firstWord, GeneData[0], articleText);
                if geneidx>0:
                    #we found a gene symbol at this word
                    #if it's a new gene for this article, add it to the list of matches
                    if not (hasSymb[geneidx] or hasName[geneidx]): 
                        matchidxList.append(geneidx);
                    hasSymb[geneidx] = True;
                    #go to next word in the article
                    articleText = articleText[len(GeneData[0][geneidx])+1:];
                    continue;
                    
                    
                #nothing was found at this word: go to next word in the article
                articleText = articleText[len(firstWord)+1:];
            
            
            #for each result, change the result from just gene index to a tuple of:
            # (article index, gene index, type of match)
            #type of match is 1=symbol, 2=name, 3=both
            for b in range(p,len(matchidxList)): 
                tm = ((1 if hasSymb[matchidxList[b]] else 0)+
                      (2 if hasName[matchidxList[b]] else 0));
                matchidxList[b] = (a,matchidxList[b],tm);
            
            #update loading bar every time we finish an article
            updateProgressBar(a+1,("Finding genes: Processed "+str(a+1)
                              +" articles out of "+str(len(ArticleData))));
        print("Finished finding genes.");
    
    trackProgress(findGenes, len(ArticleData)); #make a loading bar
    time4 = time.time();
    print("Time elapsed during gene finding: "+str(time4-time3)+" seconds.");
    
    
    ###############################################################################
    
    
    
    print("Building output...");
    
    output = {'Article':[""]*len(matchidxList),
              'Gene':[""]*len(matchidxList),
              'Gene status':[""]*len(matchidxList),
              'Type of match':[""]*len(matchidxList),
              'Locus group':[""]*len(matchidxList)};
    
    def buildOutput(updateProgressBar):
        def getStatus(b): #the status of each gene
            if needsFlag(GeneData,b): return "Shared Gene Symbol";
            if GeneData[2][b]=="APPROVED": return "Approved";
            if GeneData[2][b]=="ENTRY WITHDRAWN": return "Entry Withdrawn";
            #symbol withdrawn always says "symbol withdrawn, see ____"  
            if GeneData[2][b]=="SYMBOL WITHDRAWN": return GeneData[1][b][17:];
            return " ";
        def getType(t): #what kind of match was it?
            if t==3: return "Symbol and Name";
            if t==2: return "Symbol";
            if t==1: return "Name";
            return " ";
        for c in range(len(matchidxList)): #loop through all the matches that were found
            a = matchidxList[c][0]; #index of the article
            b = matchidxList[c][1]; #index of the gene
            output['Article'][c] = ArticleData[a];
            output['Gene'][c] = GeneData[0][b];
            output['Gene status'][c] = getStatus(b);
            output['Type of match'][c] = getType(matchidxList[c][2]);
            output['Locus group'][c] = GeneData[3][b];
            
            updateProgressBar(c+1,("Building output: Processed "+str(c+1)
                              +" outputs out of "+str(len(matchidxList))));
        print("Finished building output.");
    
    trackProgress(buildOutput,len(matchidxList)); #make a loading bar
    
    if writeAfter: 
        #save output to file
        #they entered resultsfilename in the beginning
        dfO = pd.DataFrame(data=output);
        writer = pd.ExcelWriter(resultsfilename,engine='xlsxwriter');
        dfO.to_excel(writer,sheet_name='Sheet1',index=False);
        print("Wrote results to file: "+resultsfilename);
        
    
        
    time5 = time.time();
    print("Time elapsed while building output: "+str(time5-time4)+" seconds.");

    
    
    #############################################################################
    
    
    
    def doGeneAnalysis(fileGene, fileArticle, GeneData, matchidxList, lenA,lenG):
        
        minthreshold = int(input('For charts 1,2,5,6,9,10, what should be the '+
                                 'minimum number of matches for a gene to be on the graph? '));
        numthreshold = int(input('For charts 3,4,8,9, how many of the top genes '+
                             'do you want on the graph? '));
        writeAfter = input('Do you want to save a ranked gene list to file? '+
                           '1 for yes, 0 for no: ')=='1';
        if writeAfter: 
            fileResults = input('Enter name of file to save results to: ');
        
        saveNoMatchGenes = input('Do you want to save a list of genes that have '+
                                 'no matches? 1 for yes, 0 for no: ')=='1';
        if saveNoMatchGenes: 
            fileNMG = input('Enter name of file to save genes with no matches to: ');
        
        plt.rcdefaults();
        GeneSymbols = GeneData[0]; #for bar charts
        
        # indexandsort takes a list of the number of matches as input,
        # and returns a sorted list, along with the corresponding indexes
        # example:
        # input: [0,0,1,0,0,4] = number of matches per gene
        # output:  [[5,2,4,3,1,0],  [4,1,0,0,0,0]] = [indexes of the genes, matches per gene]
        def indexandsort(countlist):
            indexes = list(range(lenG));
            indexes =  [a for _,a in sorted(zip(countlist,indexes),reverse=True)];
            countlist = sorted(countlist,reverse=True);
            return (indexes,countlist);
        
        # doubleindexandsort does the same thing, and it sorts by countlist1
        # example:
        # input: [0,0,1,3,0,4], [0,0,1,2,0,1]
        # output: [[5,3,2,4,1,0], [4,3,1,0,0,0], [1,2,1,0,0,0]]
        def doubleindexandsort(countlist1, countlist2):
            indexes = list(range(lenG));
            indexes = [a for _,a in 
                       sorted(zip(countlist1,indexes),reverse=True)];
            countlist2 = [a for _,a in 
                     sorted(zip(countlist1,countlist2),reverse=True)];
            countlist1 = sorted(countlist1,reverse=True);
            return (indexes,countlist1,countlist2);
        
        # makes a bar chart from sorted data (counts is an output of indexandsort)
        # thresh is either the number of genes to display, 
        # or the minimum number of matches for a gene to be displayed
        # which one it is depends of the True/False value threshIsNumToDisplay
        # returns the raw data for the bar chart in the form [x values, y values]
        def makeBarChart(counts, thresh, 
                       threshIsNumToDisplay, sortAlphabetical, graphtitle):
            if not threshIsNumToDisplay:
                #thresh is min hits for a gene to be displayed
                a = 0;
                while a<lenG and counts[1][a]>=thresh:
                    a = a+1;
                thresh = a;
                #now thresh is number of genes to display
                
            yvalues = counts[1][:thresh];
            xvalues = [""]*thresh;
            for a in range(thresh):
                xvalues[a] = GeneSymbols[counts[0][a]];
            
            order = list(range(thresh));
            if sortAlphabetical:
                order = [a for _,a in sorted(zip(xvalues,order))];
            
            plt.bar(order, yvalues, align='center');
            plt.xticks(order, xvalues, rotation='70'); 
                #rotation is text rotation (labels of the bar chart)
            plt.xlabel('HGNC symbol');
            plt.ylabel('Number of articles');
            plt.title(graphtitle);
            plt.show();
            return ([xvalues[a] for a in order],[yvalues[a] for a in order])
        
        # makes a grouped ber chart for sorted data (counts is an output of doubleindexandsort)
        # see comment about makeBarChart for explanation of thresh
        # returns the raw data for the bar chart in the form [x values, y1 values, y2 values]
        def makeDoubleBarChart(counts, thresh, threshIsNumToDisplay, 
                            sortAlphabetical, graphtitle, graphlegend):
            if not threshIsNumToDisplay:
                #thresh is min hits for a gene to be displayed
                a = 0;
                while a<lenG and counts[1][a]>=thresh:
                    a = a+1;
                thresh = a;
                #now thresh is number of genes to display
                
            y1values = counts[1][:thresh];
            y2values = counts[2][:thresh];
            xvalues = [""]*thresh;
            for a in range(thresh):
                xvalues[a] = GeneSymbols[counts[0][a]];
        
            order = list(range(thresh));
            if sortAlphabetical:
                order = [a for _,a in sorted(zip(xvalues,order))];
            order = np.array(order);
        
            ax = plt.subplots()[1];
            width = 0.35;
            p1 = ax.bar(order,         y1values, width, color='C0');
            p2 = ax.bar(order + width, y2values, width, color='C1');
            ax.set_xlabel('HGNC symbol');
            ax.set_ylabel('Number of articles');
            ax.set_title(graphtitle);
            ax.set_xticks(order + width / 2);
            ax.set_xticklabels(xvalues);
            ax.legend((p1[0], p2[0]), graphlegend);
            ax.autoscale_view();
            plt.show();
            return ([xvalues[a] for a in order],
                    [y1values[a] for a in order],
                    [y2values[a] for a in order]);
                    
        # makes a histogram of counts. 
        # Counts is just an array of the number of genes per article,
        # or the number of articles per gene
        # bins is just the array [0.5, 1.5, 2.5, 3.5 ....], 
        # which lets the bin/bucket of each bar correspond to an integer
        # (otherwise python makes bins of size 2.5 or so)
        def makeHistogram(counts, bins, percentFormatter, title,xlabel,ylabel):
            thehist = plt.hist(counts, bins=bins,edgecolor='black',linewidth=1);
            if percentFormatter!=None:
                plt.gca().yaxis.set_major_formatter(FuncFormatter(percentFormatter));
            plt.title(title);
            plt.xlabel(xlabel);
            plt.ylabel(ylabel);
            plt.show();
            return (list(range(1,len(thehist[0])+1)),
                    [int(a) for a in thehist[0]] if percentFormatter==None else
                    [percentFormatter(a,0) for a in thehist[0]]);
        
        
        #creation of graphs starts here
        
        countSN = [0]*lenG; #counts of all matches
        countON = [0]*lenG; #counts of matches where gene name was metioned
        countperArticle = [0]*lenA; 
        for match in matchidxList:
            countSN[match[1]] = countSN[match[1]]+1;
            if(match[2]>1):
                countON[match[1]] = countON[match[1]]+1;
            countperArticle[match[0]] = countperArticle[match[0]]+1;
        
        #sorted data - ready to make graph with
        sortedSN = indexandsort(countSN);
        sortedON = indexandsort(countON);
        sortedCombine =  doubleindexandsort(countSN,countON);
        
        #if we included 0 in the histogram, the 0s would overwhelm the actual data
        #about 99% of genes are 0s
        countsGwithout0s = [x for x in countSN if x != 0];
        countsAwithout0s = [x for x in countperArticle if x != 0];
        
        #bins make sure x-axis corresponds to the integers - 
        #otherwise python makes bins of size 2.5 or so
        binsG = [a+0.5 for a in list(range(max(countsGwithout0s)+1))];
        binsA = [a+0.5 for a in list(range(max(countsAwithout0s)+1))];
        
        #percent formatters change y axis to percent
        to_gene_percent =    lambda y,p: str(100 * y/lenG)[:7] + '%';
        to_article_percent = lambda y,p: str(100 * y/lenA)[:7] + '%';
    
        print("Gene file: "+fileGene+" ("+str(lenG)+" genes)");
        print("Article file: "+fileArticle+" ("+str(lenA)+" articles)");
        print("Number of genes mentioned in any article: "+str(len(countsGwithout0s)));
        print("Number of articles that mention any gene: "+str(len(countsAwithout0s)));
        
        # ( ,sortedSN, , , ) means all matches are counted
        # ( ,sortedON,    , , ) means only matches with gene names, or both symbol and name, are counted
        # ( ,sumCombine,      , , ) creates a double bar graph combining those two
        # ( , ,minthreshold,False, ) means it's counting all genes that have more than x hits
        # ( , ,numthreshold,True, ) means it's counting the top x genes
        # (                  ,False) means in the bar graph, genes are ordered by rank
        # (                  ,True) means in the bar graph, genes are ordered alphabetically
        fig1 = makeBarChart(sortedSN,    minthreshold,False,False,
                'Number of articles each gene symbol/name is mentioned in');
        fig2 = makeBarChart(sortedON,    minthreshold,False,False,
                'Number of articles each gene name is mentioned in');
        fig3 = makeBarChart(sortedSN,    numthreshold,True,False,
                'Number of articles each gene symbol/name is mentioned in');
        fig4 = makeBarChart(sortedON,    numthreshold,True,False,
                'Number of articles each gene name is mentioned in');
        fig5 = makeBarChart(sortedSN,    minthreshold,False,True,
                'Number of articles each gene symbol/name is mentioned in');
        fig6 = makeBarChart(sortedON,    minthreshold,False,True,
                'Number of articles each gene name is mentioned in');
        fig7 =  makeBarChart(sortedSN,    numthreshold,True,True,
                'Number of articles each gene symbol/name is mentioned in');
        fig8 = makeBarChart(sortedON,    numthreshold,True,True,
                'Number of articles each gene name is mentioned in');
        fig9 = makeDoubleBarChart(sortedCombine, minthreshold,False,False,
                'Number of articles each gene is mentioned in',
                ('Mentions by symbol or name', 'Mentions by name'));
        figA = makeDoubleBarChart(sortedCombine, minthreshold,False,True,
                'Number of articles each gene is mentioned in',
                ('Mentions by symbol or name', 'Mentions by name'));
    
        figB = makeHistogram(countsAwithout0s,binsA,None,
                      "Number of articles that mention x genes "+
                      "(articles with 0 matches excluded)",
                      'Number of genes','Number of articles');
    
        figC = makeHistogram(countsAwithout0s,binsA,to_article_percent,
                      "Percent of articles that mention x genes "+
                      "(articles with 0 matches excluded)",
                      'Number of genes','Percent of articles');
        
        figD = makeHistogram(countsGwithout0s,binsG,None,
                      "Number of genes that are mentioned in x articles "+
                      "(genes with 0 matches excluded)",
                      'Number of articles','Number of genes');
        
        figE = makeHistogram(countsGwithout0s,binsG,to_gene_percent,
                      "Percent of genes that are mentioned in x articles "+
                      "(genes with 0 matches excluded)",
                      'Number of articles','Percent of genes');
        
        #note: for the last histogram,
        #this is percentage of genes that were processed, not all genes - 
        #so if you selected locus groups, this is percentage of the genes 
        #IN THE GROUPS YOU SELECTED, not all 46000 or so
    
        
        if(writeAfter): 
            # ranks the genes by number of total matches ('total articles mentioned in')
            # ranked data was stored in sortedCombine, and was also used to make grouped bar charts
            rankedGenesList = {
                'Symbol':[""]*lenG,
                'Name':[""]*lenG,
                'Total articles mentioned in':[""]*lenG,
                'Articles that mentioned by gene name':[""]*lenG
            };
            for a in range(lenG):
                rankedGenesList['Symbol'][a] = GeneData[0][sortedCombine[0][a]];
                rankedGenesList['Name'][a]   = GeneData[1][sortedCombine[0][a]];
                rankedGenesList['Total articles mentioned in'][a] = sortedCombine[1][a];
                rankedGenesList['Articles that mentioned by gene name'][a] = sortedCombine[2][a];
                    
            dfO = pd.DataFrame(data=rankedGenesList);
            writer = pd.ExcelWriter(fileResults,engine='xlsxwriter');
            dfO.to_excel(writer,sheet_name='Sheet1',index=False);
            print("Wrote ranked genes list to file: "+fileResults);
            
        if(saveNoMatchGenes):
            #all the genes that were never mentioned at all
            nmoutput = {
                'Approved symbol':[GeneData[0][a] for a in range(lenG) if countSN[a]==0], 
                'Approved name':  [GeneData[1][a] for a in range(lenG) if countSN[a]==0], 
                'Status':         [GeneData[2][a] for a in range(lenG) if countSN[a]==0],
                'Locus group':    [GeneData[3][a] for a in range(lenG) if countSN[a]==0]
            };
            dfO = pd.DataFrame(data=nmoutput);
            writer = pd.ExcelWriter(fileNMG,engine='xlsxwriter');
            dfO.to_excel(writer,sheet_name='Sheet1',index=False);
            print("Wrote genes with no matches to file: "+fileNMG);
        
        # figures stores the raw data behind all the graphs
        # each fig object is a tuple (xvalues, yvalues)
        # except fig9 and figA, which are (xvalues,y1values,y2values)
        figures = [fig1,fig2,fig3,fig4,fig5,fig6,fig7,
                   fig8,fig9,figA,figB,figC,figD,figE];
                   
        print("Gene analysis finished!");
        return figures;
        
    # now, user can run GeneAnalysis without having to input all of the parameters
    # of doGeneAnalysis
    def GeneAnalysis():
        return doGeneAnalysis(genefilename,articlefilename,GeneData,
                            matchidxList,len(ArticleData),len(GeneData[0]));
        
      
                      
    print("...");
    print("GeneFinder finished!!!");   
    return (ArticleData,GeneData,matchidxList,output,GeneAnalysis);
    #GeneFinder also returns the GeneAnalysis function, so user can run it later
    
    
#this line begins the entire program
(ArticleData,GeneData,matchidxList,output,GeneAnalysis) = GeneFinder();
print("To run gene analysis, which creates graphs of the data, type "+
      "figures = GeneAnalysis() into the console. "+
      "This can be done multiple times.");`


    
    




```


### Support or Contact

Having trouble with Pages? Check out our [documentation](https://help.github.com/categories/github-pages-basics/) or [contact support](https://github.com/contact) and we’ll help you sort it out.
