
import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedList;


///////////// PROTEIN CLASS: CONTAINING PROTEIN NAME AND ALL ITS DIRECT NEIGHBORS   //////
public class Protein {
    String pname;
    LinkedList<Node> neighbours=new LinkedList<>();
Protein()
{}
    
Protein(String name)
{pname=name;
 //neighbours=new LinkedList<>();
} 

///////////// TO FIND THE DEGREE OF A PROTEIN ///////////////
public double pdegree()
{double degree=0.0;
 for(Node node:neighbours)
 { degree=degree+node.weight;
 }
 return degree;
}
///////////// CLONING A PROTEIN OBJECT ///////////////
public Protein proteinCopy()
{Protein newProtein=new Protein();
 newProtein.pname=this.pname;
 for(Node tnode:this.neighbours)
 { Node newNode=new Node();
   newNode.nname=tnode.nname;
   newNode.weight=tnode.weight;
   newProtein.neighbours.add(newNode);
 }
   return newProtein;
}

//////// TO FIND THE TOTAL WEIGHT OF A PROTEIN /////////////
public double findWeight()
{ double tweight=0.0;
    for(Node tempNode:neighbours)
    {tweight=tweight+tempNode.weight;
    }
 return tweight;
}
  
}

////////// COMPLEX : COLLECTION OF PROTEINS NAMES THAT FORMS A COMPLEX ////////
class Complex
{ArrayList<String> cProtein=new ArrayList<>();
 ArrayList<String>core=new ArrayList<>();
 ArrayList<String>attachment=new ArrayList<>();
}

///////////////PAIR: A PAIR REPRESENTS AN INTERACTION /////////
class Pair
{String pair1;
 String pair2;
 double weight;
}

/////////// NODE : A NODE REPRESENTS THE NEIGHBOR OF A NODE: CONTAINING THE NAME AND WEIGHT OF ITS INTERACTION ////////
class Node
{
String nname;
double weight;
}

/////////// REPRESENTS AN INTERACTION ALONG WITH THE CORRESPONDING WEIGHT///////////
class Interactions
{ 
    String protein1;
    String protein2;
    String value;
}

////////// A GENE REPRESENTS THE GENE AND ITS CORRESPONDING EXPRESSION VALUES //////////
class Gene{
String gname;
ArrayList<Double> expValue=new ArrayList<>();
boolean checkActive()
{double prev=0.0;
 double current=0.0;
 boolean active=false;
 for(double gvalue:expValue)
 { current=gvalue;
   if(prev>=0.7&&current>=0.7)
   { active=true;
     break;
   }
   else
   { prev=current;
   }
 }
 if(active)
 {return true;}
 else
 {return false;}
}

}

////////// GRAPH : TO REPRESENT A CORE CONTAINING MANY PROTEINS   ////////////////////////
class Graph {
ArrayList<Protein> ProteinChain=new ArrayList<>();

public Graph graphCopy()
{ Graph tgraph=new Graph();
  for(Protein protein:ProteinChain)
  { Protein tProtein;
    tProtein=protein.proteinCopy();
    tgraph.ProteinChain.add(protein);
  }
return tgraph;
}

////////////////// TO EVALUATE THE AVEREAGE DEGREE OF A GRAPH OR CORE  //////////
public double avgdegree()
{ double result=0.0;
  double size=ProteinChain.size();
  for(Protein protein:ProteinChain)
  { result=result+protein.pdegree();
  }
  result=result/size;
  return result;
}

////////// TO FIND THE DENSITY OF A CORE ////////////////
public double find_density()
{   double density=0.0;
    double tweight=0.0;
    double n=ProteinChain.size();
    for (Protein ProteinChain1 : ProteinChain) 
    { double pweight=0.0;
      for(Node tempNode:ProteinChain1.neighbours)
      {pweight=pweight+tempNode.weight;
      }
      tweight=tweight+pweight;

    }
    density= tweight/(n*(n-1));
 return density;
}
///////////////// ARRANGING A GRAPH AFTER REMOVAL OF VERTEX //////////////
public void arrangeGraph()
{LinkedList<String> tlist=new LinkedList<>();
 for(Protein protein:ProteinChain)
 { tlist.add(protein.pname);
 }
for(Protein protein:ProteinChain)
 { LinkedList<Node> nodeList=new LinkedList<>();
   for(Node node:protein.neighbours)
   { if(tlist.contains(node.nname))
   { nodeList.add(node);
   }
   }
   protein.neighbours=nodeList;
 }
}
////////////////// ARRANGING A GRAPH AFTER REMOVAL OF VERTEX(SIMILAR TO [arrangeGraph()])  ////////////
public Graph arrange()
{   ArrayList<String> graphProteins=new ArrayList<>();
    ArrayList<Pair> pairList=new ArrayList<>();
    Graph graph=new Graph();
    
   for(Protein protein:ProteinChain)
   { String nameP=protein.pname;
     graphProteins.add(nameP);
     for(Node node:protein.neighbours)
     { Pair tpair=new Pair();
       tpair.pair1=nameP;
       tpair.pair2=node.nname;
       tpair.weight=node.weight;
       pairList.add(tpair);
    } 
   }

     for(String nameP:graphProteins)
     { Protein protein=new Protein();
       protein.pname=nameP;
       ArrayList<String> tempNames=new ArrayList<>();
       for(Pair pair:pairList)
       { if(pair.pair1.contentEquals(nameP))
       { if(graphProteins.contains(pair.pair2))
           if(!tempNames.contains(pair.pair2))
           { Node node=new Node();
             node.nname=pair.pair2;
             node.weight=pair.weight;
             protein.neighbours.add(node);
             tempNames.add(node.nname);
            }
       }
       else if(pair.pair2.contentEquals(nameP))
       { if(graphProteins.contains(pair.pair1))
        {if(!tempNames.contains(pair.pair1))
        { Node node=new Node();
          node.nname=pair.pair1;
          node.weight=pair.weight;
          protein.neighbours.add(node);
          tempNames.add(node.nname);
        }}
       }
       }
       graph.ProteinChain.add(protein);
     }
return graph;
}



//////// removing a protein with minimum weight grom a graph ///////////
public void removeMinWeight()
{Protein minWeightProtein=ProteinChain.get(0);
 for(Protein protein:ProteinChain)
 { if(protein.findWeight()<minWeightProtein.findWeight())
   { minWeightProtein=protein;
   }
 }
 ProteinChain.remove(minWeightProtein);
}
}
