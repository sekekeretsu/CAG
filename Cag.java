
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.charset.StandardCharsets;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.Scanner;
import java.util.StringTokenizer;

/**
 *
 * @author seke14
 */
public class Cag {
      static String weighted_data;//="C:\\\\coach\\\\gavin.txt";
      static String gene_data;//="C:\\\\coach\\\\gep\\\\gene.txt";
      static String benchmark_data;//="C:\\\\coach\\\\benchmark\\\\sgd.txt";
     static  double density_th;
     static double similarity_th;
   
    public static void main(String[] args) throws IOException {
          
    weighted_data=args[0];
     gene_data=args[1];
     benchmark_data=args[2];
     density_th=Double.parseDouble(args[3]);
     similarity_th=Double.parseDouble(args[4]);
      ProteinComplex PC=new ProteinComplex();
       PC.distinctProteins(weighted_data);
       PC.readBenchmarkComplex(benchmark_data);
       PC.gepreader(gene_data);
       PC.matchProtienToGene();
       PC.removeInvalidInteraction();
       PC.writeInteractions();
   //  PC.checkRedundentInteraction();                       
   //  PC.displayGenes();
       PC.proteinNeighborChain();
   //  PC.displayProteinNeighbour();
       PC.similarNeighbours(similarity_th);
   //  PC.displaysimilarNeighbour();
   //  PC.testProteinNeighbour();
       PC.coreComplex(density_th);
   
       PC.redundancyFilter();
   //  PC.coreDisplay();
   //  PC.addingAttachmentByGene();
       PC.addingAttachments();
   //  PC.displayPredictedComplex();
  
       PC.CompareToFindDetectedComplex();
       PC.CompareToFindRealComplex();
  //   PC.displayComplex(); 
       PC.writeDetectedComplex();
       PC.displayComplexCollectionDetails();
       PC.evaluationDisplay();
    }
}
class ProteinComplex{
   ArrayList<Interactions> InteractionList=new ArrayList<>();
   ArrayList<Interactions> InteractionList_filtered=new ArrayList<>();
   ArrayList<String> dist_element=new ArrayList<>();
   ArrayList<Protein> neighbourList=new ArrayList<>();
   ArrayList<Gene>   geneList=new ArrayList<>();
   ArrayList<String> dist_genename=new ArrayList<>();
   ArrayList<String> dist_protein=new ArrayList<>();
   ArrayList<String> nomatch=new ArrayList<>();
   ArrayList<Graph> similarNeighbourList=new ArrayList<>();
   ArrayList<Graph> coreList=new ArrayList<>();
   ArrayList<Graph> coreListFiltered=new ArrayList<>();
   ArrayList<Complex> predictedComplexList=new ArrayList<>();
   ArrayList<Complex> benchmarkComplexList=new ArrayList<>();
   ArrayList<Complex> matchComplexList=new ArrayList<>(); 
   ArrayList<Complex> matchRealComplexList=new ArrayList<>();
   ArrayList<Complex> coveredRealComplex=new ArrayList<>();
   ArrayList<Complex> coveredRealComplex1=new ArrayList<>();
   ArrayList<Complex> coreComplexList=new ArrayList<>();
        double detected;
        double predicted;
        double benchmark;
        double realMatched;
        double pre;
        double recall;
        double fmeasure;
        int no_protein;
        int no_interaction;
        int no_gene;
      
        
    ///////////////////////////   INPUT: READING INPUT PPI NETWORK DATA    /////////////////////////////////
    
    public  void distinctProteins(String aFileName) throws IOException {
       
         fFilePath=Paths.get(aFileName);
         try (Scanner scanner =  new Scanner(fFilePath, ENCODING.name())){
         while (scanner.hasNextLine())
         {        Interactions tInteraction=new Interactions();
                
             String aLine=scanner.nextLine();
                         String name1=new String();
                         String name2=new String();
                         String value=new String();
                         StringTokenizer st=new StringTokenizer(aLine);
                     if (st.hasMoreTokens())
                         name1 = st.nextToken();
                      if(st.hasMoreTokens())
                         name2 = st.nextToken();
                      if(st.hasMoreTokens())
                         value=st.nextToken();
                      tInteraction.protein1=name1;
                      tInteraction.protein2=name2;
                      tInteraction.value=value;
                      InteractionList.add(tInteraction);
                  if(!dist_element.contains(name1.trim()))
                      dist_element.add(name1.trim());
                  if(!dist_element.contains(name2.trim()))
                      dist_element.add(name2.trim());
         
         }System.out.print(" end line ");
       }  
}
     
        
   /////////////////////////INPUT : READING REFERENCE OR BENCHMARK DATA  //////////////////////////////////////////
   
    public void readBenchmarkComplex(String bname) throws IOException
    { fFilePath_benchmark=Paths.get(bname);
     try(Scanner scanner=new Scanner(fFilePath_benchmark,ENCODING.name()))
     { 
     while(scanner.hasNextLine())
     { String aLine;aLine=scanner.nextLine();
       Complex tComplex=new Complex();   
     StringTokenizer stzer=new StringTokenizer(aLine);
        while(stzer.hasMoreTokens())
        {String thisString=stzer.nextToken();
         tComplex.cProtein.add(thisString);
        } benchmarkComplexList.add(tComplex);
     }
     }
    }
      //////////////// INPUT DATA : READING GENE EXPRESSION PROFILE INPUT DATA /////////////
     public void gepreader(String fFileName) throws IOException
    {
       fFilePath_GeneExp=Paths.get(fFileName);
       try (Scanner scanner =  new Scanner(fFilePath_GeneExp, ENCODING.name())){
         while (scanner.hasNextLine()){
                Gene tempGene=new Gene();
                 String aLine=scanner.nextLine();
                 StringTokenizer st = new StringTokenizer(aLine);
                 String firstToken;
                 firstToken=st.nextToken();
                 String nameToken=st.nextToken();
                 dist_genename.add(nameToken);
                 tempGene.gname=nameToken;
                   while (st.hasMoreTokens()) {
                   tempGene.expValue.add(Double.parseDouble(st.nextToken()));
           }  
                  geneList.add(tempGene);
      }  
    }
      
    }
      
       ///////////  FINDING PROTEIN HAVING SIMILAR EXPRESION FROM A PROTEIN NEIGHBORHOOD ////////////////
 
  public void similarNeighbours(double similarity_th)
   {
       for(int x=0;x<neighbourList.size();x++)
       { Graph tGraph=new Graph();
         Protein protein;//=new Protein();
         protein=neighbourList.get(x);
         protein=protein.proteinCopy();    //.......proteinCopy().../////
         String pname=protein.pname;
         Gene gene1;//=new Gene();
         gene1=findGene(pname);
         if(gene1==null)
         {System.out.println(" gene exp not found @similarneighbour gene1");
         continue;
         }
         for(int y=0;y<protein.neighbours.size();y++)
         {  String pname2=protein.neighbours.get(y).nname;
            Gene gene2=findGene(pname2);//new Gene();
            if(gene2==null)
            {   System.out.println(" gene exp not found @@similarneighbour gene2");
                continue;}
            else 
            { Double similarity=evalSimilarity(gene1,gene2);
              if(similarity>similarity_th)
              {   Protein protein2=findProtein(pname2);
                 if(protein2==null)
                 { System.out.println("     this protein in similar Neighbour cannot be empty");
                 }
                 else
                 {
                  tGraph.ProteinChain.add(protein2);
                 }
                 }
            }
         }
         tGraph.ProteinChain.add(protein);
         similarNeighbourList.add(tGraph);
         }
         System.out.println("size of the similar neighbour list: "+ similarNeighbourList.size());
           }
 
      //////////////// FINDING THE DIRECT NEGHBORS OF A PROTEIN   /////////////////////////////////////
 
   LinkedList<Node> findNeighbors(String passedProteinname) throws IOException
    {  LinkedList<Node> templl=new LinkedList<>();
       LinkedList<String> tempNeighbours=new LinkedList<>();
        for(int x=0;x<InteractionList_filtered.size();x++)
        { String name1=InteractionList_filtered.get(x).protein1;
          String name2=InteractionList_filtered.get(x).protein2;
          String value=InteractionList_filtered.get(x).value;
           
          if(name1.contentEquals(name2))
          { continue;
          }
         
          if(passedProteinname.contentEquals(name1))
          {
              if(tempNeighbours.contains(name2))
                {continue;}
              else    
                 {   Node tnode=new Node();
                     tnode.nname=name2;
                     tnode.weight=Double.parseDouble(value);
                     templl.add(tnode);
                     tempNeighbours.add(name2);
                 }
          }
          else if(passedProteinname.contentEquals(name2))
          {
              if(tempNeighbours.contains(name1))
              {continue;}
              else
              {
              Node tnode=new Node();
              tnode.nname=name1;
              tnode.weight=Double.parseDouble(value);
              templl.add(tnode);
              tempNeighbours.add(name1);
             }
          }
        }
        return templl;
    }
    
   
   ///////////////  IDENTIFYING COMPLEX CORES:   //////////////////////////////////
   public void coreComplex(double density_thr)
   {  
       for(int x=0;x<similarNeighbourList.size();x++)
      { Graph tempGraph=similarNeighbourList.get(x);
       tempGraph.arrangeGraph();
    //  tempGraph=tempGraph.graphPruning();
       if(tempGraph.ProteinChain.size()>=2)
       {  
           double d=tempGraph.find_density();
           if(d>=density_thr)
           { coreList.add(tempGraph);
           }
   else
           { LinkedList<Graph> subGraphList=new LinkedList<>();
            subGraphList=removeCore(tempGraph,density_thr);
            for(Graph tGraph:subGraphList)
            { 
            if(tGraph.find_density()>density_thr)
            {   tGraph=tGraph.arrange();
                coreList.add(tGraph);
            }
           else
           { 
           while(d<density_thr&&(tGraph.ProteinChain.size()>2))
           { tGraph.removeMinWeight();
             tGraph=tGraph.arrange();
             d=tGraph.find_density();
           }
           if(tGraph.ProteinChain.size()>=2&&d>=density_thr)
           {   tGraph=tGraph.arrange();
               coreList.add(tGraph);
           }}
        }
     }
   }}
   }
   
 /////////////// SPLIT AND MERGE LOW DENSITY CORES TO GET SUB-COMPONENTS OR SUB CORES WITH HIGH DENSITY ////////////////////  
 LinkedList<Graph> removeCore(Graph graph,double den)
 {LinkedList<Graph> result=new LinkedList<>();
  if(graph.find_density()>den) 
  { result.add(graph);
  }
  else if(graph.ProteinChain.size()<=2)
  {
     result.add(graph);
  }
  else
  { LinkedList<Graph> graphList=new LinkedList<>();
    LinkedList<Protein> Lower=new LinkedList<>();
    LinkedList<Protein> Higher=new LinkedList<>();
    LinkedList<String> hString=new LinkedList<>();
    double avgDegree=graph.avgdegree();
            for(Protein protein:graph.ProteinChain)
            {
            if(protein.pdegree()>=avgDegree)
            { Higher.add(protein);
              hString.add(protein.pname);
            }
            else
           {Lower.add(protein);
            }
            }
            
    for(Protein protein:Lower)
    { LinkedList<Node> nodeList=new LinkedList<>();
    for(Node node:protein.neighbours)
    { if(hString.contains(node.nname))
    {}
    else
    {
    nodeList.add(node);
    }
    }protein.neighbours=nodeList;
    }
    for(Protein protein:Lower)
    { int flag=0;
       LinkedList<String> node_names=new LinkedList<>();
       node_names.add(protein.pname);
      for(Node tnode:protein.neighbours)
      { node_names.add(tnode.nname);
      }
      for(Graph tgraph:graphList)
      {LinkedList<String> graphPNames=new LinkedList<>();
      for(Protein tprotein:tgraph.ProteinChain)
      { graphPNames.add(tprotein.pname);
        for(Node tnode:tprotein.neighbours)
        { graphPNames.add(tnode.nname);
        }
      }
      for(String name1:node_names)
      { if(graphPNames.contains(name1))
      { flag=1;
      tgraph.ProteinChain.add(protein);
      break;
      }
      }
      if(flag==1)
      { break;
      }
      }
      if(flag==0)
      { Graph newGraph=new Graph();
        newGraph.ProteinChain.add(protein);
        graphList.add(newGraph);
      }
    }
    LinkedList<Graph> preResult=new LinkedList<>();
    for(Graph tgraph:graphList)
    { LinkedList<Graph> tempResult=new LinkedList<>();
      tempResult=removeCore(tgraph,den);
      preResult.addAll(tempResult);
    }
    for(Graph tGraph:preResult)
    { tGraph.ProteinChain.addAll(Higher);
      tGraph.arrangeGraph();
      result.add(tGraph);
    }   
  }
 return result;
 }
 
  ///////////////////////////////  ADDING ATTACHMENT PROTEINS TO THE COMPLEX CORE  //////////////////////////////////////////
  
public void addingAttachments()
  {
   for(Graph tgraph:coreListFiltered)
   {Graph graph1=new Graph();
    Complex newComplex=new Complex();
    LinkedList<String> coreProteinList=new LinkedList<>();
    LinkedList<String> attachment=new LinkedList<>();
    double coreListSize=tgraph.ProteinChain.size();
    for(Protein protein1:tgraph.ProteinChain)
    { coreProteinList.add(protein1.pname);
      newComplex.cProtein.add(protein1.pname);
      Protein protein2=findProtein(protein1.pname);
      protein2=protein2.proteinCopy();  
      graph1.ProteinChain.add(protein2);
    }
    LinkedList<String> candidateProteins=new LinkedList<>();
    for(Protein protein1:graph1.ProteinChain)
    { for(Node node1:protein1.neighbours)
     { if(candidateProteins.contains(node1.nname)||coreProteinList.contains(node1.nname))
        { continue;
        }
     else
     { candidateProteins.add(node1.nname);
     }
     }
    }
    for(String pname:candidateProteins)
    { double intersection=0.0;
      for(Protein protein:graph1.ProteinChain)
      { for(Node node1:protein.neighbours)
        { if(pname.contentEquals(node1.nname))
        { intersection++;
          break;
        }
        }
      }
      if((intersection*intersection)>=(0.5*coreListSize))
      {attachment.add(pname);
       newComplex.cProtein.add(pname);
      }
    }
    predictedComplexList.add(newComplex);
   }
  }
  
   /////////////////////////////////////////  FILTERING REDUNDANT CORES //////////////////////////////////////
 
  public void redundancyFilter()
   { coreListFiltered.add(coreList.get(0));
      for(Graph tcore1:coreList)
      { double NA_max=0.0;
        Graph max=new Graph();
        for(Graph tcore2:coreListFiltered)
        { double na;
        na=neighbourAffinity(tcore1,tcore2);
        if(na>NA_max)
        {
            NA_max=na;
            max=tcore2;
        }
        else
        {continue;
        }
        
        }
        if(NA_max<0.8)
        { coreListFiltered.add(tcore1);
        }
        else
        {
            double den1=tcore1.find_density();
            double den2=max.find_density();
            double size1=tcore1.ProteinChain.size();
            double size2=max.ProteinChain.size();
            if((den1*size1)>=(den2*size2))
            { coreListFiltered.add(tcore1);
              coreListFiltered.remove(max);
            }
           }
          }
       System.out.println("size of filtered list:"+coreListFiltered.size());
   }
 
  ////////////////  COMPUTING THE NEIGHBORHOOD AFFINITY OF TWO COMPLEX CORES TOBE USED TO REMOVE REDUNDANCY //////////////
  
   double neighbourAffinity(Graph graph1,Graph graph2)
   {double result=0.0;
   double size1=graph1.ProteinChain.size();
   double size2=graph2.ProteinChain.size();
   double match=0.0;
   if(size1>=size2)
   {  LinkedList<String> ProteinList=new LinkedList<>();
      for(Protein protein1:graph1.ProteinChain)
      {ProteinList.add(protein1.pname);
      }
   for(Protein protein2:graph2.ProteinChain)
      { if(ProteinList.contains(protein2.pname))
      { match++;
      }
      }
   }
   else
   {
       LinkedList<String> ProteinList2=new LinkedList<>();
       for(Protein protein2:graph2.ProteinChain)
       {
           ProteinList2.add(protein2.pname);
       }
       for(Protein protein1:graph1.ProteinChain)
       { if(ProteinList2.contains(protein1.pname))
       { match++;
       }
       }
   }
   double denom=size1*size2;
   double NA=(match*match)/denom;
   result=NA;
    return result;
   }
   
 
 
 
     //////////////////////// ....PREPROCESS INPUT : TO CHECK REDUNDANT INTERACTIONS ... ///////////////////////   
  
  public void checkRedundentInteraction()
    { int count=0;
    for(Interactions interaction:InteractionList_filtered)
    {String name1=interaction.protein1;
     String name2=interaction.protein2;
     for(Interactions interaction2:InteractionList_filtered)
     {if(interaction2.protein1.contentEquals(name2)&&interaction2.protein2.contentEquals(name1))
     { System.out.println("there is redundency:"+(count++)+" "+name1+" "+" "+name2);break;}}
    }
    }
    
    /////////////////...... WRITING INTERACTIONS OF PROTEINS HAVING EXPRESSION VALUE ... //////////////////////////
 
   public void writeInteractions()
    { try{
    File f=new File("output_filteredPPI.txt");
    f.createNewFile();
    FileWriter Fw=new FileWriter(f);
    for(Interactions interaction:InteractionList_filtered)
    {Fw.write(interaction.protein1+"\t"+interaction.protein2+"\t"+"1");
     Fw.write("\n");
    }
    Fw.flush();
    Fw.close();
    }
   catch(Exception e){
            System.out.println("EXCEPTION"+e.getMessage());
        }
    System.out.println(" the size of the old [[after]] filtering Interaction LIst is:"+ InteractionList.size());
   System.out.println(" the size of the new Interaction LIst is                    :"+ InteractionList_filtered.size());
  
   System.out.println("total proteins in the network :"+dist_element.size());
   System.out.println("total proteins in the network witout a matching gene :"+nomatch.size());
   System.out.println("total proteins in the network with a match:"+dist_protein.size());
    }
    
    /////////////////////////////////////......... TO DISPLAY PREDICTED COMPLEXES .....///////////////////////////
    public void displayPredictedComplex()
    { System.out.println(".................displaying predicted complex:...................  "); 
        for(Complex tComplex:predictedComplexList)
    { for(String tprotein:tComplex.cProtein)
      {System.out.print(tprotein+"\t"); 
      }
    System.out.println();
    } 
    }
    
   
    
    ///////////////////......  DISPLAYS  PREDICTED,DETECTED ,REFERENCE COMPLEX ETC. ...... //////////////////
      void displayComplexCollectionDetails()
     { System.out.println();
       System.out.println("real complex size                 :"+benchmarkComplexList.size());
       System.out.println("predicted complex size            :"+predictedComplexList.size());
       System.out.println("matching complex size             :"+matchComplexList.size());
       System.out.println("real complex that has match match :"+coveredRealComplex.size());
       System.out.println("size of unfiltered core list size :"+coreList.size());
       System.out.println("size of filtered core list size   :"+coreListFiltered.size());
    
     }
         public void writeDetectedComplex() throws IOException
       {
       File f=new File("output_DetectedComplex.txt");
       f.createNewFile();
       FileWriter fw=new FileWriter(f);
       
       System.out.println();
         System.out.println();
           System.out.println();
       System.out.println("oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo");
       System.out.println("ooooooooooooooooooooooo  WRITING AND PRINTING DETECTED COMPLEXES TO FILE  oooooooooooooooooooooooooooooooooooooooooo");
       System.out.println("oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo");
         System.out.println();
           System.out.println();
       for(Complex complex:matchComplexList)
      { for(String protein:complex.cProtein)
        { System.out.print(protein);
          fw.write(protein);
        System.out.print("\t");
        fw.write("\t");}
      System.out.println();
      fw.write("\n");
      }
       fw.flush();
       fw.close();
        
   }
   /////////////////////////   DISPLAYING PREDICTD COMPLEXES THAT HAS A MATCH IN THE REFERENCE DATA   //////////////
   public void displayComplex()
   {
       
      System.out.println();
         System.out.println();
           System.out.println();
       System.out.println("oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo");
       System.out.println("ooooooooooooooooooooooooooooo  PRINTING DETECTED COMPLEXES TO FILE  oooooooooooooooooooooooooooooooooooooooooooooooo");
       System.out.println("oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo");
         System.out.println();
           System.out.println();
           for(Complex complex:matchComplexList)
      { for(String protein:complex.cProtein)
        { System.out.print(protein);
        System.out.print("\t");}
      System.out.println();
      }
       
   }
   
   
  
////////////////////////// TO REMOVE INTERACTIONS CONTAINING PROTEINS WITH NO EXPRESSION VALUE  ///////////////////////////////////
   public void removeInvalidInteraction()
   {   int count=0;
       System.out.println(" the size of the old [[before]] filtering Interaction LIst is:"+ InteractionList.size());
       for(int x=0;x<InteractionList.size();x++)
        {
            if(nomatch.contains(InteractionList.get(x).protein1)||nomatch.contains(InteractionList.get(x).protein2))
                {count++;continue;}
            else
                {if(InteractionList.get(x).protein1.contentEquals(InteractionList.get(x).protein2))
                {continue;}
                else
                {InteractionList_filtered.add(InteractionList.get(x));
                }}
        }
   System.out.println(" the size of the old [[after]] filtering Interaction LIst is:"+ InteractionList.size());
   System.out.println(" the size of the new Interaction LIst is                    :"+ InteractionList_filtered.size());
   System.out.println(" count of invalid interaction                               : "+count);
   System.out.println("total proteins in the network :"+dist_element.size());
   System.out.println("total proteins in the network witout a matching gene :"+nomatch.size());
   System.out.println("total proteins in the network with a match:"+dist_protein.size());
   for(Interactions i1:InteractionList_filtered)
   {  
       String name1=i1.protein1;
       String name2=i1.protein2;
      
          if(!(dist_protein.contains(name1)))
          dist_protein.add(name1);
          if(!(dist_protein.contains(name2)))
          dist_protein.add(name2);
   }
   }

   
   ////////////////////////    DIPLAYING PROTEIN WITH HIGH SIMILARITY VALUE IN A SIMILAR NEIGHBORHOOD ///////////////////////
   public void displaysimilarNeighbour()
   {
   System.out.println("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX");
   System.out.println("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX");
   System.out.println("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX");
   System.out.println("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX");
   System.out.println("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX");
  for(int x=0;x<similarNeighbourList.size();x++)
  { Graph tempgraph=similarNeighbourList.get(x);
     System.out.println(" "+ tempgraph.ProteinChain.size());
      for(int y=0;y<tempgraph.ProteinChain.size();y++)
  {   Protein tempProtein=tempgraph.ProteinChain.get(y);
          if(tempProtein!=null)
          System.out.print(" "+tempProtein.pname);
                  else 
              System.out.println("accessing a null  ");
  }   System.out.println();
  }
   }
   
    
   ///////// FINDING PROTEIN : RETURNS A PROTIEN NAME AND ITS DIRECTLY CONNECTED NEIGHBORS   //////////////
   Protein findProtein(String pname)
   {Protein protein=new Protein();
    for(int x=0;x<neighbourList.size();x++)
     {   
        if(neighbourList.get(x).pname.contentEquals(pname))
        {
         protein= neighbourList.get(x);
         break;
        }
        if(x==(neighbourList.size()-1)&&protein.pname==null)
        {System.out.println("the protein is not found:"+ pname);
           }
        }
    Protein proteinR;
    proteinR=protein.proteinCopy();
    return proteinR;
    }
   
   ////////////// ............. FINDS A GENE HAVING EXPRESSION VALUES......//////////////
   Gene findGene(String gName)
   { Gene returnGene=null;
     for(int x=0;x<geneList.size();x++)
     { Gene tgene=geneList.get(x);
       if(tgene.gname.contentEquals(gName))
       { returnGene=tgene;break;} //////////// break???
     }
   return returnGene;
   }
   
   ///////////// EVALUATES THE SIMILARITY IN THE EXPRESSION OF TWO PROTEINS : USING PEARSONS CORRELATION COEFFICIENT //////////////////////
  Double evalSimilarity(Gene gene1,Gene gene2)
   {Double result;
    int size=36;
    Double x,y;
    Double xSum=0.0;
    Double ySum=0.0;
    Double xSumSquare=0.0;
    Double ySumSquare=0.0;
    Double num,Denom,Denom1;
    Double xy=0.0;
    for(int i=0;i<size;i++)
    {
        x=gene1.expValue.get(i);
        y=gene2.expValue.get(i);
        xSum=xSum+x;
        ySum=ySum+y;
        xSumSquare=xSumSquare+(x*x);
        ySumSquare=ySumSquare+(y*y);
        xy=xy+(x*y);
    }
    num=(size*xy)-(xSum*ySum);
    Denom1=((size*xSumSquare)-(xSum*xSum))*((size*(ySumSquare))-(ySum*ySum));
    Denom=Math.sqrt(Denom1);
    result=num/Denom;
    
   // System.out.println("num:"+num);
    
    //System.out.println("denom:"+Denom);
    System.out.println("PCC : "+result);
    return result;
   }
   
  //////////// FINDING PROTEINS THAT HAS CORRESPONDING GENE EXPRESSION PROFILE ////////////////////
   public void matchProtienToGene()
   { for(int x=0;x<dist_element.size();x++)
   { String temp=dist_element.get(x);
     if(!(dist_genename.contains(temp)))
      { nomatch.add(temp);
        }
   }
   }
   
   //////////////////////////// TO DISPLAY GENE AND CORRESPONDING EXPRESSION VALUES //////////////
   
     void displayGenes()
     { for(Gene a:geneList)
     { 
        System.out.print(a.gname);
        for(Double b:a.expValue)
        {System.out.print(" "+b);}
                System.out.println();
     }
          System.out.println("the sizeof geneList or total genes from GEP       : "+geneList.size());
          System.out.println("the sizeof distinct elements from weightd PPI     : "+dist_element.size());
          System.out.println("the sizeof  dist gene from GEP                    : "+dist_genename.size());
          System.out.println("the sizeof proteins form PPI mathing genes from GEP: "+dist_protein.size());
          System.out.println("the sizeof Proteins with no match or expression   : "+nomatch.size());
     }
    
    
     //////////////////////////////// TO DISPLAY PROTEIN-PROTEIN INTERACTIONS /////////////////////////////////////
     void displayIntractions()
     {   System.out.println("size of the interactions:"+InteractionList.size());
         System.out.println("size of the interactions:"+dist_element.size());
         for(int i=0;i<InteractionList.size();i++)
         { 
          System.out.println(i+" "+InteractionList.get(i).protein1+" "+InteractionList.get(i).protein2+" "+InteractionList.get(i).value);
         }
     }
     
     //////////////////   Creates a LINKEDLIST of PROTEIN ; These Proteins contain a protein name AND Protein neighbors /////////////////////
     //////////////////   and other nodes that intereact with the Protein                         /////////////////////
    public void proteinNeighborChain() throws IOException
    {
     for(int i=0;i<dist_protein.size();i++)
    {
         String tdist_protein=dist_protein.get(i);
         Protein protein=new Protein(tdist_protein);
         protein.neighbours=findNeighbors(tdist_protein);    
         neighbourList.add(protein);
    
    }
    }
   
    ///////////////////////////// TO DIPLAY THE INDIVIDUAL AND ITS CORRESPONDING NEIGHBOR OR PROTEIN THAT HAS INTERACTION WITH IT///////////////
    void displayProteinNeighbour()
    {  for(int i=0;i<neighbourList.size();i++)
    {System.out.println( "PROTEIN :"+neighbourList.get(i).pname);
         for(int j=0;j<neighbourList.get(i).neighbours.size();j++)
         {System.out.print( neighbourList.get(i).neighbours.get(j).nname+" ");}
    System.out.println( "\n");
    }
    }
    
  
    ///////////////////////////////// FINDING REAL COMPLEXES THAT HAS A MATCH WITH ONE OR MORE PREDICTED COMPLEXES  ////////////////////////////////////////
     public void CompareToFindRealComplex()
    { System.out.println("........ FINDING REAL COMPLEXES THAT HAS MATCH IN THE PREDICTED COMPLEX SET   ..........");
  
        for(Complex pComplex:benchmarkComplexList)
        { LinkedList<String> pComplexElements=new LinkedList<>();
           double size1=pComplex.cProtein.size();
                for(String string:pComplex.cProtein)
                 { pComplexElements.add(string);
                 }
             
        Complex maxComplex=new Complex();
        double maxCloseness=0.0;
    for(Complex bComplex:predictedComplexList)
        {    double match=0.0;
             for(String string1:bComplex.cProtein)
              {    
                 if(pComplexElements.contains(string1))
                 {     match++;
                 }
             }
             double size2=bComplex.cProtein.size();
             double prod1=match*match;
             double prod2=size1*size2;
             double closeness=prod1/prod2;
             if(closeness>maxCloseness)
             { maxCloseness=closeness;
               maxComplex=bComplex;
             }
        }
    if(maxCloseness>0.255)
    { matchRealComplexList.add(pComplex);
      coveredRealComplex1.add(maxComplex);
    }
    }
    }
     
     //////////////////////////// DETECTED COMPLEXES : FINDING PREDICTED COMPLEXES THAT HAS ATLEAST A MATCH IN THE REFERENCE COMPLEX SET     ///////////////////////////////
   
    public void CompareToFindDetectedComplex()
    { System.out.println(".......  DETECTED COMPLEXES:  FINDING PREDICTED COMPLEXES THAT HAS A MATCH IN THE REFERENCE COMPLEX SET   .......");
        for(Complex pComplex:predictedComplexList)
           { LinkedList<String> pComplexElements=new LinkedList<>();
             double size1=pComplex.cProtein.size();
             for(String string:pComplex.cProtein)
                 { pComplexElements.add(string);
                 }
             
        Complex maxComplex=new Complex();
        double maxCloseness=0.0;
    for(Complex bComplex:benchmarkComplexList)
        {    double match=0.0;
             for(String string1:bComplex.cProtein)
              {    
                 if(pComplexElements.contains(string1))
                 {     match++;
                 }
             }
             double size2=bComplex.cProtein.size();
             double prod1=match*match;
             double prod2=size1*size2;
             double closeness=prod1/prod2;
             if(closeness>maxCloseness)
             { maxCloseness=closeness;
               maxComplex=bComplex;
             }
        }
    if(maxCloseness>0.255)
    { matchComplexList.add(pComplex);
      coveredRealComplex.add(maxComplex);
    }
    }
    }
    
    //////////////////////////////// TO DISPLAY THE COMPLEX CORES PREDICTED ///////////////////////////////  
    public void coreDisplay()
   { System.out.println();
     for(Graph tgraph:coreListFiltered)
     { for(Protein tprotein:tgraph.ProteinChain)
        {System.out.print(tprotein.pname+":\t");
          for(Node node:tprotein.neighbours)
          {System.out.print(node.nname+ " ");
          }System.out.println();
          
        }System.out.println();
   }
   }
    
     ///////////////...... EVALUATE   AND DISPLAY FINAL EVALUEATED VALUES OF PRECISION REACALL AND FMEASURE ...../////////////
    public void evaluationDisplay()
    {
         detected=0;
         predicted=0;
         benchmark=0;
         realMatched=0;
         pre=0;
         recall=0;
        benchmark=benchmarkComplexList.size();
        predicted=predictedComplexList.size();
        detected=matchComplexList.size();
        realMatched=matchRealComplexList.size();
        benchmark=benchmarkComplexList.size();
        pre=detected/predicted;
        recall=realMatched/benchmark;
        fmeasure=2*(pre*recall)/(pre+recall);
        no_protein=dist_protein.size();
        no_interaction=InteractionList_filtered.size();
        no_gene=dist_genename.size();
        System.out.println();
        System.out.println("............... FINAL EVALUATED VALUE OF PRECISION RECALL AND F-MEASURE................");
        System.out.println();
        System.out.println("total proteins in the network having expression data : "+dist_protein.size());
        System.out.println("total no. of the Protein-Protein Interactions        : "+ InteractionList_filtered.size());

        System.out.println();
        System.out.println("no of real complexes in reference dataset    : "+benchmarkComplexList.size());
        System.out.println("no of total predicted complex                : "+predictedComplexList.size());
        System.out.println("no of predicted complexes that matched a real complex in reference Dataset    : "+matchComplexList.size());
        System.out.println("no of real complexes having a match in the Predicted complex set              : "+matchRealComplexList.size());
        System.out.println("                precision :"+pre);
        System.out.println("                    recall:"+recall);
        System.out.println("                f-measure :"+fmeasure);
        System.out.println();
        System.out.println();
        
    }
  

    private Path fFilePath_benchmark;
    private Path fFilePath_GeneExp;
    private static Path fFilePath;
    private final static Charset ENCODING = StandardCharsets.UTF_8;
}////..............CLASS CLOSED HERE.............//////
