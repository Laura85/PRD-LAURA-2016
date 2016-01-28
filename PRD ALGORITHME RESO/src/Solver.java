import ilog.concert.IloException;
import ilog.concert.IloLinearNumExpr;
import ilog.concert.IloNumVar;
import ilog.cplex.IloCplex;


public class Solver {
	public void solveProblem(){
		int n = 4;
		int m = 3;
		int V = 2;
		int T = 10;
		int cV = 1000;
		
		//int[][] matrice=new int[n][m];
		
		int[][] p =
		    {
		        { 5, 3, 9 } , // tableau [0] de int
		        { 4, 2, 5 } ,
		        { 2, 6, 1 } ,
		        { 1, 5, 2 } 
		    };
		int[][] x =
		    {
		        { 0,1,0,0 } , // tableau [0] de int
		        { 0,0,0,0 } ,
		        { 0,0,0,1 } ,
		        { 1,0,0,0 } 
		    };
		int[][] y =
		    {
		        { 0,0,1,1 } , // tableau [0] de int
		        { 1,1,0,0 } 
		    };
		
		int[] d ={ 25,26,16,19 };
		int[] hWIP ={ 2,4,1,3 };
		int[] hFIN ={ 2,1,4,3 };
		int[] piM ={ 10,20,30,40 };
		int[] q ={ 50,100,200,100 };
		int[] F ={ 30,60 };
		
		
		try {
			IloCplex cplex = new IloCplex();
			
			//variables
			IloNumVar[][] C = new IloNumVar[m][];
			for (int i=0;i<m;i++){
				C[i] = cplex.numVarArray(n, 0, Double.MAX_VALUE);
			}
			IloNumVar[] PT = cplex.numVarArray(n, 0, Double.MAX_VALUE);
			IloNumVar[] ICFin = cplex.numVarArray(n, Double.MIN_VALUE, Double.MAX_VALUE);
			IloNumVar[] ICWip = cplex.numVarArray(n, Double.MIN_VALUE, Double.MAX_VALUE);
			IloNumVar IC = cplex.numVar(Double.MIN_VALUE, Double.MIN_VALUE);
			
			//objectif
			IloLinearNumExpr obj = cplex.linearNumExpr();
			obj.addTerm(1.0, IC);
			//obj.addTerm(1.0, V*cV);
			//TODO ajouter V cV à l'objectif
			for(int j=0;j<n;j++){
				obj.addTerm(piM[j],PT[j]);
			}
			cplex.addMinimize(obj);
			
			//contraintes
			
			//contr_ress_1
			for (int i=0;i<m;i++){
				for (int j1=0;j1<n;j1++){
					for (int j2=0;j2<n;j2++){
						if(j1!=j2){
							cplex.addGe(C[i][j2], cplex.prod(cplex.sum(C[i][j1],p[i][j2]), x[j1][j2]));
						}
					}
				}
			}
			
			//contr_fin
			for (int j=0;j<n;j++){
				int cst = 0;
				for(int v=0;v<V;v++){
					cst += (F[v]*y[j][v]);
				}
				cplex.addLe(C[m][j], cst);
			}
			
			//contr_fin2
			for (int j=0;j<n;j++){
				cplex.addGe(C[1][j], p[1][j]);
			}
			
			//contr_gamme
			for (int j=0;j<n;j++){
				for (int i=0;i<m-1;i++){
					cplex.addGe(C[i+1][j], cplex.sum(C[i][j],p[i+1][j]));
				}
			}
			
			//calcul_ICFIN
			for (int j=0;j<n;j++){
				int cst = 0;
				for(int v=0;v<V;v++){
					cst += (F[v]*y[j][v]);
				}
				cplex.addEq(ICFin[j], cplex.prod(cplex.prod(cplex.sum(cst,cplex.prod(-1.0, C[m][j])),q[j]), hFIN[j]));
			}
			
			//calcul_ICWIP
			for (int j=0;j<n;j++){
				cplex.addEq(ICWip[j], cplex.prod(cplex.prod(cplex.sum(C[m][j],cplex.prod(-1.0, C[1][j])),q[j]), hWIP[j]));
			}
			
			//calcul_ICTOT
			IloLinearNumExpr expr = cplex.linearNumExpr();
			for(int j=0;j<n;j++){
				expr.addTerm(1.0,ICFin[j]);
				expr.addTerm(1.0,ICWip[j]);
			}
			
			//calcul_retard
			for (int j=0;j<n;j++){
				int cst = 0;
				for(int v=0;v<V;v++){
					cst += (F[v]*y[j][v]);
				}
				cplex.addGe(PT[j], cst + T - d[j]);
			}
			
			//solve problem
			cplex.solve();
			//end
		} catch (IloException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
}
