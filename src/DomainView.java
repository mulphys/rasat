/*
	Domain Viewer
	Author: Andrei Smirnov (andrei.v.smirnov@gmail.com)
*/

/*
	To get the Render mode switch uncomment all lines with comboBox
*/

import java.util.*;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import javax.swing.*;

import java.applet.Applet;
import java.awt.Image;
import java.awt.Event;
import java.awt.Graphics;
import java.awt.Dimension;
import java.io.*;
import java.net.URL;
import java.util.Hashtable;
import java.awt.image.IndexColorModel;
import java.awt.image.ColorModel;
import java.awt.image.MemoryImageSource;
import java.awt.event.*;
import java.awt.*;
import java.util.Vector;

///class	DomainObject
///{	public String	name=null;
///	Color	color;
///	public int	hits=0;
///}	

class Statistics
{
	float
		probability=0.0F,
		error=1.0F;
	public void setProbability(float p){ probability = p; }
	public float getProbability(){ return probability; }
	public void setError(float e){ error = e; }
	public float getError(){ return error; }
}

class ObjData
{	final static int DIM=3;
	static int nsample=0,
		idata, //index of the data set
		Dim[]=new int[DIM],
		X0[]=new int[DIM],
		dX[]=new int[DIM],
		mangles=0,//max and actual angles
		hits[][][][],
		angles[][][][];
	static float dangle=1.0F;//default angle
	static void setDataIndex(int i){ idata=i; }
	static int getDataIndex(){ return idata; }
	static public int getSampleSize(){ return nsample; }
	static public int[] getDim(){ return Dim; }
	static public int getDim(int i){ return Dim[i]; }
	static public int[] getX0() { return X0; }
	static public int getX0(int i) { return X0[i]; }
	static public int[] getAngles(int i, int j, int k)
	{	return angles[i][j][k];
	}
	static public int[] getHits(int i, int j, int k)
	{	return hits[i][j][k];
	}
	static void getHitStatistics
	(	int DIM, 
		double dY[], //size of one cell
		double DY[], //size of the domain
		double y[], //normalized position of the point (0...1)
		double angle,
		Statistics statistics
	) // RETURNS HIT-STATISTICS
	{	int IX[]=new int[DIM];
		for(int i=0;i<DIM;i++)
		{	//Cell-index of the point
			IX[i]=(int)Math.floor(DY[i]*y[i]/dY[i]-X0[i]);
			if(IX[i]<0||IX[i]>=Dim[i])
			{	statistics.setProbability(-1.0F);
				return;//data not available
			}
		}
		int
			a[]=getAngles(IX[0],IX[1],IX[2]),
			nangles=a[0],
			h[]=getHits(IX[0],IX[1],IX[2]);
		//Binary search
		int k0=0,k1=nangles;
		while(k1-k0>1)
		{	int k=(k0+k1)/2;
			if(angle>a[k]*dangle)
				k0=k;
			else
				k1=k;
		}
		double p=
		(double)h[k0];
		if(k0<k1)
		{	//Linear interpolation:
			double da=dangle*(a[k1]-a[k0]);
			p+=
			+ ((double)(h[k1]-h[k0])*(angle-a[k0]*dangle))
			/ da;
		}
		statistics.setProbability((float)(p/(double)nsample));
		statistics.setError((float)(2.0/Math.sqrt(p)));
	}
//-	static float getHitProbability
//-	(	int DIM, 
//-		double dY[], //size of one cell
//-		double DY[], //size of the domain
//-		double y[], //normalized position of the point (0...1)
//-		double angle
//-	) // RETURNS HIT-PROBABILITY
//-	{	int IX[]=new int[DIM];
//-		for(int i=0;i<DIM;i++)
//-		{	//Cell-index of the point
//-			IX[i]=(int)Math.floor(DY[i]*y[i]/dY[i]-X0[i]);
//-			if(IX[i]<0||IX[i]>=Dim[i])
//-				return -1.0F;//data not available
//-		}
//-		int
//-			a[]=getAngles(IX[0],IX[1],IX[2]),
//-			nangles=a[0],
//-			h[]=getHits(IX[0],IX[1],IX[2]);
//-		//Binary search
//-		int k0=0,k1=nangles;
//-		while(k1-k0>1)
//-		{	int k=(k0+k1)/2;
//-			if(angle>a[k]*dangle)
//-				k0=k;
//-			else
//-				k1=k;
//-		}
//-		double p=
//-		(double)h[k0];
//-		if(k0<k1)
//-		{	//Linear interpolation:
//-			double da=dangle*(a[k1]-a[k0]);
//-			p+=
//-			+ ((double)(h[k1]-h[k0])*(angle-a[k0]*dangle))
//-			/ da;
//-		}
//-		return (float)(p/(double)nsample);
//-	}
//+	static float getHitProbabilityInterp
//+	(	int DIM, 
//+		double dY[], //size of one cell
//+		double DY[], //size of the domain
//+		double y[], //normalized position of the point (0...1)
//+		double angle
//+	) // RETURNS HIT-PROBABILITY INTERPOLATING FROM THE COARSE GRID
//+	{	int IX[]=new int[DIM];
//+		float 
//+			Y[]=new float[DIM],
//+			D[]=new float[DIM];
//+		for(int i=0;i<DIM;i++)
//+		{	D[i]=(float)(dX[i]*dY[i]); //size of the coarse-grid cell
//+			Y[i]=(float)(DY[i]*y[i]);//absolute position of the point
//+		}
//+		for(int i=0;i<DIM;i++) 
//+		{	IX[i]=(int)((Math.floor(Y[i]/dY[i])-X0[i])/dX[i]);
//+			if(IX[i]<0||IX[i]>=Dim[i])
//+				return -1.0F;//data not available
//+		}
//+		//RETRIEVE HITS FROM SAMPLE DATA:
//+		int ni[]=new int[DIM];
//+		for(int i=0;i<DIM;i++)
//+		{	ni[i]=2;
//+			if(IX[i]+1>Dim[i]) ni[i]=1;
//+		}
//+		float	P[][][]=new float[2][2][2];//Hit probabilities
//+		for(int i0=0;i0<ni[0];i0++)
//+		{	int j0=IX[0]+i0;
//+		for(int i1=0;i1<ni[1];i1++)
//+		{	int j1=IX[1]+i1;
//+		for(int i2=0;i2<ni[2];i2++)
//+		{	int j2=IX[2]+i2,
//+				a[]=getAngles(j0,j1,j2),
//+				nangles=a[0],
//+				h[]=getHits(j0,j1,j2);
//+			//Binary search
//+			int k0=0,k1=nangles;
//+			while(k1-k0>1)
//+			{	int k=(k0+k1)/2;
//+				if(angle>a[k]*dangle)
//+					k0=k;
//+				else
//+					k1=k;
//+			}
//+			float da=dangle*(float)(a[k1]-a[k0]);
//+			//INTERPOLATE FOR ANGLE:
//+			P[i0][i1][i2]=
//+			(float)h[k0]
//+			+ (float)((h[k1]-h[k0])*(angle-a[k0]*dangle))
//+			/ (float)da;
//+		}
//+		}
//+		}
//+		//INTERPOLATE FOR XYZ:
//+		float d[]=new float[DIM];//differential quotes
//+		for(int i=0;i<DIM;i++)d[i]=(float)(Y[i]-dY[i]*(IX[i]*dX[i]+X0[i]))/D[i];
//+		float p=0.0F;
//+		float p2[][]=new float[2][2];
//+		for(int i0=0;i0<ni[0];i0++)
//+		for(int i1=0;i1<ni[1];i1++)
//+		{	p2[i0][i1]=P[i0][i1][0];
//+			if(ni[2]==2)
//+				p2[i0][i1]+=(P[i0][i1][1]-P[i0][i1][0])*d[0];
//+		}
//+		float p1[]=new float[2];
//+		for (int i0=0;i0<ni[0];i0++)
//+		{	p1[i0]=p2[i0][0];
//+			if(ni[1]==2)p1[i0]+=(p2[i0][1]-p2[i0][0])*d[1];
//+		}
//+		p=p1[0];
//+		if(ni[0]==2)
//+			p+=(p1[1]-p1[0])*d[1];
//+		return p/(float)nsample;
//+	}
	static public void Load(int iobj, URL url)
	{	try
		{	String filename="obj"+iobj+".dat.gz"; //object data file name
			System.out.println("Loading sample data from: "+filename);
			StreamTokenizer st = new StreamTokenizer
			(	new BufferedReader
				(	new InputStreamReader
					(	new GZIPInputStream
						(///	new FileInputStream(filename)
							new URL(url,filename).openStream()
						)
					)
				)
			);
//			st.eolIsSignificant(true);
//			st.commentChar('#');
			int ix[] = new int[DIM];
			int 
				iangle[] = null,
				anglehits[]=null;
			while(st.nextToken() != StreamTokenizer.TT_EOF)
			{	switch (st.ttype) 
				{	default:
					break;
					case StreamTokenizer.TT_EOL:
					break;
					case StreamTokenizer.TT_WORD:
					if ("SAMPLE".equals(st.sval)) 
					{	st.nextToken();
						if (st.nextToken() != StreamTokenizer.TT_NUMBER) 
						throw new Exception("Format exception: "+st.toString());
						nsample=(int) st.nval;	
					}
					else
					if (st.sval.equals("NX"))
					{	st.nextToken();
						for(int i=0;i<DIM;i++)
						{	if (st.nextToken() != StreamTokenizer.TT_NUMBER) 
							throw new Exception("Format exception: "+st.toString());
							Dim[i]=(int) st.nval;	
						}
						hits=new int[Dim[0]][Dim[1]][Dim[2]][];
						angles=new int[Dim[0]][Dim[1]][Dim[2]][];
					}
					else
					if (st.sval.equals("X0"))
					{	st.nextToken();
						for(int i=0;i<DIM;i++)
						{	if (st.nextToken() != StreamTokenizer.TT_NUMBER) 
							throw new Exception("Format exception: "+st.toString());
							X0[i]=(int) st.nval;	
						}
					}
					else
					if (st.sval.equals("DX"))
					{	st.nextToken();
						for(int i=0;i<DIM;i++)
						{	if (st.nextToken() != StreamTokenizer.TT_NUMBER) 
							throw new Exception("Format exception: "+st.toString());
							dX[i]=(int) st.nval;	
						}
					}
					else
					if (st.sval.equals("NANGLES"))
					{	st.nextToken();
						if (st.nextToken() != StreamTokenizer.TT_NUMBER) 
						throw new Exception("Format exception: "+st.toString());
						mangles=(int) st.nval;
						iangle = new int[mangles];
						anglehits = new int[mangles];
					}
					else
					if (st.sval.equals("DANGLE"))
					{	st.nextToken();
						if (st.nextToken() != StreamTokenizer.TT_NUMBER) 
						throw new Exception("Format exception: "+st.toString());
						dangle=(float) st.nval;
					}
					else
					if (st.sval.equals("XYZ"))
					{	st.nextToken();
						for(int i=0;i<DIM;i++)
						{	st.nextToken();
							ix[i]=(int)st.nval;
						}
						int nangles=0;
///						st.nextToken();
						while(st.nextToken()==StreamTokenizer.TT_NUMBER)
						{	iangle[nangles]=(int)st.nval;
							if (st.nextToken() != StreamTokenizer.TT_NUMBER) 
							throw new Exception("Format exception: "+st.toString());
							anglehits[nangles]=(int)st.nval;
							if(nangles++>=mangles)
							{	System.out.println
								("Too many angles: "+nangles+" >= "+mangles);
								break;
							}
						}
						hits[ix[0]][ix[1]][ix[2]]=new int[nangles+1];
						hits[ix[0]][ix[1]][ix[2]][0]=nangles;
						angles[ix[0]][ix[1]][ix[2]]=new int[nangles+1];
						angles[ix[0]][ix[1]][ix[2]][0]=nangles;
						for(int i=0;i<nangles;i++)
						{	angles[ix[0]][ix[1]][ix[2]][i+1]=iangle[i];
							hits[ix[0]][ix[1]][ix[2]][i+1]=anglehits[i];
						}
					}
					break;
///					else 
///					{	while 
///						(	st.nextToken() != StreamTokenizer.TT_EOL
///							&& st.ttype != StreamTokenizer.TT_EOF
///						);
///					}
				}// end switch
///				if (st.ttype != StreamTokenizer.TT_EOF)
///				throw new Exception(st.toString());
			}
///					st.close();
			//CONSISTENCY CHECK:
			for(int i=0;i<Dim[0];i++)
			for(int j=0;j<Dim[1];j++)
			for(int k=0;k<Dim[2];k++)
			{	if(angles[i][j][k]==null||hits[i][j][k]==null)
				throw new Exception("CONSISTENCY CHECK FAILED AT GRID NODE: "+i+","+j+","+k);
			}
		}// end try
		catch (IOException e)
		{	System.out.println("IO-Exception: "+e);
		}
		catch (Exception  e)
		{	System.out.println("File format exception: "+e);
		}
	}
	static public int Load(int iobj, int[] IX, URL url)
	{	if(IX==null) return -1;
		String indexFileName="index.dat";
		int idata=-1;
		try //Read Index File
		{	StreamTokenizer in = new StreamTokenizer
			(	new BufferedReader
				(	new InputStreamReader
					(	new URL(url,indexFileName).openStream()
					)
				)
			);
			//Read in the components
			int[]
				Iorigin = new int[DIM],
				Ndomains = new int[DIM],
				DomDim = new int[DIM];	
			in.nextToken(); //NAME
			for(int i=0;i<DIM;i++)
			{	in.nextToken(); 
				Iorigin[i] = (int)in.nval;
			}
			in.nextToken(); //NAME
			for(int i=0;i<DIM;i++)
			{	in.nextToken(); 
				Ndomains[i] = (int)in.nval;
			}
			in.nextToken(); //NAME
			for(int i=0;i<DIM;i++)
			{	in.nextToken(); 
				DomDim[i] = (int)in.nval;
			}
			int n=1;
			for(int i=0;i<DIM;i++) n*=Ndomains[i];
			idata=0;
			for(int i=0;i<DIM;i++)
			{	int j=IX[i]>=Iorigin[i]?(IX[i]-Iorigin[i])/DomDim[i]:-1;
				if(j<0 || j>=Ndomains[i])
					throw new Exception
					("Data not available for index: "+i+", value: "+IX[i]);
				n/=Ndomains[i];
				idata+=j*n;
			}
		}
		catch(Exception e)
		{	System.out.println("Warning: Reading "+indexFileName+" file:"+ e);
			return -1;
		}
		if(idata<0) return idata;
		try //Read Data File
		{	String filename="s"+idata+"o"+iobj+".dat.gz";
			System.out.println("Loading data-set: "+filename);
			StreamTokenizer st = new StreamTokenizer
			(	new BufferedReader
				(	new InputStreamReader
					(	new GZIPInputStream
						(///	new FileInputStream(filename)
							new URL(url,filename).openStream()
						)
					)
				)
			);
//			st.eolIsSignificant(true);
//			st.commentChar('#');
			int ix[] = new int[DIM];
			int 
				iangle[] = null,
				anglehits[]=null;
			while(st.nextToken() != StreamTokenizer.TT_EOF)
			{	switch (st.ttype) 
				{	default:
					break;
					case StreamTokenizer.TT_EOL:
					break;
					case StreamTokenizer.TT_WORD:
					if ("SAMPLE".equals(st.sval)) 
					{	st.nextToken();
						if (st.nextToken() != StreamTokenizer.TT_NUMBER) 
						throw new Exception("Format exception: "+st.toString());
						nsample=(int) st.nval;	
					}
					else
					if (st.sval.equals("NX"))
					{	st.nextToken();
						for(int i=0;i<DIM;i++)
						{	if (st.nextToken() != StreamTokenizer.TT_NUMBER) 
							throw new Exception("Format exception: "+st.toString());
							Dim[i]=(int) st.nval;	
						}
						hits=new int[Dim[0]][Dim[1]][Dim[2]][];
						angles=new int[Dim[0]][Dim[1]][Dim[2]][];
					}
					else
					if (st.sval.equals("X0"))
					{	st.nextToken();
						for(int i=0;i<DIM;i++)
						{	if (st.nextToken() != StreamTokenizer.TT_NUMBER) 
							throw new Exception("Format exception: "+st.toString());
							X0[i]=(int) st.nval;	
						}
					}
					else
					if (st.sval.equals("DX"))
					{	st.nextToken();
						for(int i=0;i<DIM;i++)
						{	if (st.nextToken() != StreamTokenizer.TT_NUMBER) 
							throw new Exception("Format exception: "+st.toString());
							dX[i]=(int) st.nval;	
						}
					}
					else
					if (st.sval.equals("NANGLES"))
					{	st.nextToken();
						if (st.nextToken() != StreamTokenizer.TT_NUMBER) 
						throw new Exception("Format exception: "+st.toString());
						mangles=(int) st.nval;
						iangle = new int[mangles];
						anglehits = new int[mangles];
					}
					else
					if (st.sval.equals("DANGLE"))
					{	st.nextToken();
						if (st.nextToken() != StreamTokenizer.TT_NUMBER) 
						throw new Exception("Format exception: "+st.toString());
						dangle=(float) st.nval;
					}
					else
					if (st.sval.equals("XYZ"))
					{	st.nextToken();
						for(int i=0;i<DIM;i++)
						{	st.nextToken();
							ix[i]=(int)st.nval;
						}
						int nangles=0;
///						st.nextToken();
						while(st.nextToken()==StreamTokenizer.TT_NUMBER)
						{	iangle[nangles]=(int)st.nval;
							if (st.nextToken() != StreamTokenizer.TT_NUMBER) 
							throw new Exception("Format exception: "+st.toString());
							anglehits[nangles]=(int)st.nval;
							if(nangles++>=mangles)
							{	System.out.println
								("Too many angles: "+nangles+" >= "+mangles);
								break;
							}
						}
						hits[ix[0]][ix[1]][ix[2]]=new int[nangles+1];
						hits[ix[0]][ix[1]][ix[2]][0]=nangles;
						angles[ix[0]][ix[1]][ix[2]]=new int[nangles+1];
						angles[ix[0]][ix[1]][ix[2]][0]=nangles;
						for(int i=0;i<nangles;i++)
						{	angles[ix[0]][ix[1]][ix[2]][i+1]=iangle[i];
							hits[ix[0]][ix[1]][ix[2]][i+1]=anglehits[i];
						}
					}
					break;
///					else 
///					{	while 
///						(	st.nextToken() != StreamTokenizer.TT_EOL
///							&& st.ttype != StreamTokenizer.TT_EOF
///						);
///					}
				}// end switch
///				if (st.ttype != StreamTokenizer.TT_EOF)
///				throw new Exception(st.toString());
			}
///					st.close();
			//CONSISTENCY CHECK:
			for(int i=0;i<Dim[0];i++)
			for(int j=0;j<Dim[1];j++)
			for(int k=0;k<Dim[2];k++)
			{	if(angles[i][j][k]==null||hits[i][j][k]==null)
				throw new Exception("CONSISTENCY CHECK FAILED AT GRID NODE: "+i+","+j+","+k);
			}
			setDataIndex(idata);
		}// end try
		catch (IOException e)
		{	System.out.println("IO-Exception: "+e);
			return -1;
		}
		catch (Exception  e)
		{	System.out.println("File format exception: "+e);
			return -1;
		}
		return idata;
	}
}

/** Wireframe */

class WireFrame 
{
	static int DIM=3;
	int	dim,Dim[];
	static Color col[][];
	//Particles Injection data:
	static int ninj,minj,//number of injection points
		inj0,inj1,//injection point and the first particle 
		iframe0;//start index of the shooting frame
	static int trapped=0; //number of particels trapped
	static double[] 
		X0,dX,DX,//Domain origin and dimensions
		Xinj=null, //injection coordinates
		Vair=null; //Air velocity
	static float vairscale=3.0F;
///		Vel; //particle velocities
	static double	vdir;// velocity direction
	static int[]	Cinj; //injection points color
	static byte[] grid; //domain geometry
	//Points:
	int nvert, mvert, maxvert;
	float vert[];
	int tvert[];
///	int ZsortMap[];
	int ncon, maxcon;
	int con[];
	byte com[];
	int	nshades;
	boolean transformed;
	Matrix3D mat;
	float xmin, xmax, ymin, ymax, zmin, zmax;

///	private final static int bgGrey = 192;
///	private static int maxr;
///	private int Rl;
///	private int Gl;
///	private int Bl;
///	static 
///	{	data = new byte[R * 2 * R * 2];
///		int mr = 0;
///		for (int Y = 2 * R; --Y >= 0;) 
///		{	int x0 = (int) (Math.sqrt(R * R - (Y - R) * (Y - R)) + 0.5);
///			int p = Y * (R * 2) + R - x0;
///			for (int X = -x0; X < x0; X++) 
///			{	int x = X + hx;
///				int y = Y - R + hy;
///				int r = (int) (Math.sqrt(x * x + y * y) + 0.5);
///				if (r > mr)
///				mr = r;
///				data[p++] = r <= 0 ? 1 : (byte) r;
///			}
///		}
///		maxr = mr;
///	}
///	private final int blend(int fg, int bg, float fgfactor) 
///	{
///		return (int) (bg + (fg - bg) * fgfactor);
///	}

	WireFrame () 
	{
		mat = new Matrix3D ();
		mat.xrot(20);
		mat.yrot(30);
		ncon=0;
		nshades=16;
	}
	WireFrame (Space space, Domain domain, Vector components) 
	{	this();
		dim=space.getDim();
		Dim=domain.getDim();
		int ncomp=components.size();
		col = new Color[nshades][ncomp];
		dX = space.getdX();
		DX = domain.getDX();
		X0 = domain.getX0();
		grid=domain.getgrid();
		float[]
			dY = new float[dim],
			Y0 = new float[dim];
		int nv=0;//vertex counter
		vert=null;
		tvert=null;
///		ZsortMap=null;
		nvert=nv;
		//Set components
		for (int icomp=1;icomp<ncomp;icomp++)
		{	Components	comp=(Components)components.elementAt(icomp);
			for (int ishade = 0; ishade < nshades; ishade++) 
			{	double shade=1.0*(1-Math.pow(ishade/((float)nshades-1.0), 2.3));
				int
					grey=170,
					r=(int)((1-shade)*comp.color.getRed()+shade*grey),
					g=(int)((1-shade)*comp.color.getGreen()+shade*grey),
					b=(int)((1-shade)*comp.color.getBlue()+shade*grey);
				col[ishade][icomp] = new Color(r,g,b);
			}
			float[]
				x1 = new float[dim],
				y1 = new float[dim];
			if(comp.boundaryContours==null)continue;
			//ADD BOUNDARY CONTOURS
			for(int idir=0;idir<dim;idir++)
			{	if(comp.boundaryContours[idir]==null)continue;
				BoundaryContour	contour=comp.boundaryContours[idir];
				while (contour!=null)
				{
					x1[0]=(float)contour.getVertexIndX(0);
					x1[1]=(float)contour.getVertexIndY(0);
					x1[2]=(float)contour.getPlaneInd();
					for (int i=0;i<dim;i++)
					{	int	j=(i+dim-idir+2)%dim;
						y1[i]=(float)X0[j]+(float)dX[j]*x1[j];
					}
					addVert(y1[0], y1[1], y1[2]);nv++;
					for(int ic=1;ic<contour.Length();ic++)
					{	x1[0]=(float)contour.getVertexIndX(ic);
						x1[1]=(float)contour.getVertexIndY(ic);
						for (int i=0;i<dim;i++)
						{	int	j=(i+dim-idir+2)%dim;
							y1[i]=(float)X0[j]+(float)dX[j]*x1[j];
						}
						addVert(y1[0], y1[1], y1[2]);
						add(nv-1,nv,icomp);nv++;
					}
					contour=contour.next;
				}
			}
		}
		//Set particles color (strored as 0-component color)
		for (int ishade = 0; ishade < nshades; ishade++) 
		{	double shade=1.0*(1-Math.pow(ishade/((float)nshades-1.0), 2.3));
			int
				grey=170,
				r=(int)((1-shade)*Cinj[0]+shade*grey),
				g=(int)((1-shade)*Cinj[1]+shade*grey),
				b=(int)((1-shade)*Cinj[2]+shade*grey);
			col[ishade][0] = new Color(r,g,b);
		}
		{	//ADD NODES FOR THE SHOOTING FRAME:
			//Injection coordinates:
			float[] xinj = new float[dim];
			for(int i=0;i<dim;i++)
				xinj[i]=(float)Math.floor(DX[i]*Xinj[i]);
			iframe0=nvert;
			for(int i0=0;i0<2;i0++)
			for(int i1=0;i1<2;i1++)
			for(int i2=0;i2<2;i2++)
			{	addVert //Use the injection coordinates 
				(	xinj[0], 
					xinj[1], 
					xinj[2]
				);
				nv++;
			}
			// ADD FRAME EDGES
			for(int i=0;i<2;i++)
			for(int j=0;j<2;j++)
			{	add(iframe0+2*j,iframe0+2*j+1,0);
				add(iframe0+2*j,iframe0+4*i+2*j,0);
				add(iframe0+2*j+1,iframe0+4*i+2*j+1,0);
				add(iframe0+4*i+2*j,iframe0+4*i+2*j+1,0);
			}
			//ADD PARTICLES
			///for(int iinj=0;iinj<ninj;iinj++)
			for(int iinj=0;iinj<ninj+1;iinj++)
			{	addVert((float)xinj[0], (float)xinj[1], (float)xinj[2]);
				nv++;
			}
			inj0=nvert-ninj-1; 
			inj1=inj0+1; 
			int  jnj0=DIM*inj0,jnj1=DIM*inj1;
			for (int i=0;i<DIM;i++)
			{	vert[jnj1+i]=vert[jnj0+i]+vairscale*(float)Vair[i];
			}
			add(inj0,inj1,0);
		}
	}
	void setFrame()
	{	//Display the shooting area:
		// SET SHOOTING FRAME NODES
		for(int i=0;i<DIM;i++)
		{	float
				dx=(float)dX[i],
				Dx=ObjData.getDim(i)*dx,
				x0=ObjData.getX0(i)*dx;
			int iv=0;
			for(int i0=0;i0<2;i0++)
			for(int i1=0;i1<2;i1++)
			for(int i2=0;i2<2;i2++)
			{	int j=i==0?i0:(i==1?i1:i2);
				vert[DIM*(iframe0+iv++)+i]=x0+(float)j*Dx;
			}
		}
	}
	void hideFrame()
	{	// Hide the shooting area:
		for(int i=0;i<DIM;i++)
		{	int iv=0;
			for(int i0=0;i0<2;i0++)
			for(int i1=0;i1<2;i1++)
			for(int i2=0;i2<2;i2++)
				vert[DIM*(iframe0+iv++)+i]=(float)getinj(i);
		}
	}
	/** Create a 3D model by parsing an input stream */
///	WireFrame (InputStream is) throws IOException, FileFormatException 
///	{	this();
///		StreamTokenizer st = new StreamTokenizer(new BufferedReader(new InputStreamReader(is)));
///		st.eolIsSignificant(true);
///		st.commentChar('#');
///	scan:
///	while (true) {
///		switch (st.nextToken()) {
///			default:
///		break scan;
///			case StreamTokenizer.TT_EOL:
///		break;
///			case StreamTokenizer.TT_WORD:
///		if ("v".equals(st.sval)) {
///			double x = 0, y = 0, z = 0;
///			if (st.nextToken() == StreamTokenizer.TT_NUMBER) {
///			x = st.nval;
///			if (st.nextToken() == StreamTokenizer.TT_NUMBER) {
///				y = st.nval;
///				if (st.nextToken() == StreamTokenizer.TT_NUMBER)
///				z = st.nval;
///			}
///			}
///			addVert((float) x, (float) y, (float) z);
///			while (st.ttype != StreamTokenizer.TT_EOL &&
///				st.ttype != StreamTokenizer.TT_EOF)
///			st.nextToken();
///		} else if ("f".equals(st.sval) || "fo".equals(st.sval) || "l".equals(st.sval)) {
///			int start = -1;
///			int prev = -1;
///			int n = -1;
///			while (true)
///			if (st.nextToken() == StreamTokenizer.TT_NUMBER) {
///				n = (int) st.nval;
///				if (prev >= 0)
///				add(prev - 1, n - 1);
///				if (start < 0)
///				start = n;
///				prev = n;
///			} else if (st.ttype == '/')
///				st.nextToken();
///			else
///				break;
///			if (start >= 0)
///			add(start - 1, prev - 1);
///			if (st.ttype != StreamTokenizer.TT_EOL)
///			break scan;
///		} else {
///			while (st.nextToken() != StreamTokenizer.TT_EOL
///				&& st.ttype != StreamTokenizer.TT_EOF);
///		}
///		}
///	}
///	is.close();
///	if (st.ttype != StreamTokenizer.TT_EOF)
///		throw new FileFormatException(st.toString());
///	}
	/** Delete a vertex from this model */
	int delVert(int ivert) 
	{
		if (ivert >= nvert) return nvert;
		nvert--;
		for (int i=ivert;i<nvert;i++)
		{	int j = 3*i, j1=j+3;
			for(int k=0;k<3;k++)
				vert[j+k] = vert[j1+k];
		}
		return nvert;
	}
	/** Add a vertex to this model */
	int addVert(float x, float y, float z) 
	{
		int i = nvert;
		if (i >= maxvert)
			if (vert == null) 
			{
				maxvert = 100;
				vert = new float[maxvert * 3];
			} else 
			{	maxvert *= 2;
				float nv[] = new float[maxvert * 3];
				System.arraycopy(vert, 0, nv, 0, vert.length);
				vert = nv;
			}
		i *= 3;
		vert[i] = x;
		vert[i + 1] = y;
		vert[i + 2] = z;
		if(mvert<=nvert)mvert=nvert+1;
		return nvert++;
	}
	/** Add a line from vertex p1 to vertex p2 */
	void add(int p1, int p2, int icomp) 
	{
		int i = ncon;
		if (p1 >= nvert || p2 >= nvert)
			return;
		if (i >= maxcon)
		if (con == null) 
		{
			maxcon = 100;
			con = new int[maxcon];
			com = new byte[maxcon];
		}	else 
		{	maxcon *= 2;
			int nv[] = new int[maxcon];
			System.arraycopy(con, 0, nv, 0, con.length);
			con = nv;
			byte nc[] = new byte[maxcon];
			System.arraycopy(com, 0, nc, 0, com.length);
			com = nc;
		}
		if (p1 > p2) 
		{	int t = p1;
			p1 = p2;
			p2 = t;
		}
		con[i] = (p1 << 16) | p2;
		com[i] = (byte)icomp;
		ncon = i + 1;
	}
	/** Transform all the points in this model */
	void transform() 
	{
		if (transformed || nvert <= 0)
			return;
		if (tvert == null || tvert.length < nvert * 3)
			tvert = new int[nvert*3];
		mat.transform(vert, tvert, nvert);
		transformed = true;
	}
	float getDepth(int T)
	{
		int
			p1=((T >> 16) & 0xFFFF),
			p2=(T & 0xFFFF),
			iz1=3*p1+2,
			iz2=3*p2+2;
		float
			z1=vert[iz1],
			z2=vert[iz2],
			z=z1>z2?z2:z1;
		return -z;
	}
	/* Quick Sort implementation
	*/
	private void quickSort(int a[], byte b[], int left, int right)
	{
		int leftIndex = left;
		int rightIndex = right;
	///	int partionElement;
		float partionElement;
		if ( right > left)
		{	/* Arbitrarily establishing partition element as the midpoint of
				* the array.
				*/
			partionElement = getDepth(a[ ( left + right ) / 2 ]);
			// loop through the array until indices cross
			while( leftIndex <= rightIndex )
			{	/* find the first element that is greater than or equal to
				 * the partionElement starting from the leftIndex.
				 */
				while( ( leftIndex < right ) && ( getDepth(a[leftIndex]) < partionElement ) )
					++leftIndex;
				/* find an element that is smaller than or equal to
				 * the partionElement starting from the rightIndex.
				 */
				while( ( rightIndex > left ) &&
					( getDepth(a[rightIndex]) > partionElement ) )
					--rightIndex;
				// if the indexes have not crossed, swap
				if( leftIndex <= rightIndex )
				{	swap(a, leftIndex, rightIndex);
					swapb(b, leftIndex, rightIndex);
					++leftIndex;
					--rightIndex;
				}
			}
			/* If the right index has not reached the left side of array
				* must now sort the left partition.
				*/
			if( left < rightIndex )
				quickSort( a, b, left, rightIndex );
			/* If the left index has not reached the right side of array
				* must now sort the right partition.
				*/
			if( leftIndex < right )
				quickSort( a, b, leftIndex, right );
		}
	}
	private void swap(int a[], int i, int j)
	{
		int T;
		T = a[i];
		a[i] = a[j];
		a[j] = T;
	}
	private void swapb(byte a[], int i, int j)
	{
		byte T;
		T = a[i];
		a[i] = a[j];
		a[j] = T;
	}
	/** eliminate duplicate lines */
	void compress() 
	{
		int limit = ncon;
		int c[] = con;
		byte cm[] = com;
		quickSort(con, com, 0, ncon - 1);
		int d = 0;
		int pp1 = -1;
		for (int i = 0; i < limit; i++) 
		{	int p1 = c[i];
			byte m1 = cm[i];
			if (pp1 != p1) 
			{	c[d] = p1;
				cm[d]= m1;
				d++;
			}
			pp1 = p1;
		}
		ncon = d;
	}
	/** Paint this model to a graphics context.	It uses the matrix associated
	with this model to map from model space to screen space.
	The next version of the browser should have double buffering,
	which will make this *much* nicer */
	void paint(Graphics g) 
	{
		if (vert == null || nvert <= 0)
			return;
		transform();
		int v[] = tvert;
///		int zs[] = ZsortMap;
///		if (zs == null) 
///		{
///			ZsortMap = zs = new int[nvert];
///			for (int i = nvert; --i >= 0;)
///			zs[i] = i * 3;
///		}
///		/*
///		* I use a bubble sort since from one iteration to the next, the sort
///		* order is pretty stable, so I just use what I had last time as a
///		* "guess" of the sorted order.  With luck, this reduces O(N log N)
///		* to O(N)
///		*/
///		for (int i = nvert - 1; --i >= 0;) 
///		{	boolean flipped = false;
///			for (int j = 0; j <= i; j++) 
///			{	int a = zs[j];
///				int b = zs[j + 1];
///				if (v[a + 2] > v[b + 2]) 
///				{
///					zs[j + 1] = a;
///					zs[j] = b;
///					flipped = true;
///				}
///			}
///			if (!flipped)
///			break;
///		}
		int lg = 0;
		int lim = ncon;
		int c[] = con;
		if (lim <= 0 || nvert <= 0) return;
		for (int i = 0; i < lim; i++) 
		{	int T = c[i];
			int 
				p1=((T >> 16) & 0xFFFF) * 3,
				p2=(T & 0xFFFF) * 3,
				grey = v[p1 + 2] + v[p2 + 2];
			if (grey < 0)grey = 0;
			if (grey >= nshades)	grey = nshades-1;
			if (i*nshades+grey != lg) 
			{	lg = i*nshades+grey;
				g.setColor(col[grey][(int)com[i]]);
			}
			g.drawLine
			(	v[p1], v[p1 + 1],
				v[p2], v[p2 + 1]
			);
		}
		//Display particles:
		{	//Injection point
			////int i=nvert-ninj,
			int i=inj0,
				p=3*i,
				grey = v[p + 2];
			if (grey < 0)grey = 0;
			if (grey >= nshades)	grey = nshades-1;
			if (i*nshades+grey != lg) 
			{	lg = i*nshades+grey;
				g.setColor(col[grey][0]);
			}
			g.fillOval
			(	v[p]-4, v[p + 1]-4,
				8, 8 // SIZE OF THE POINT SOURCE
			);

		}
		//Other points
		for (int i = nvert-ninj+2; i < nvert; i++) 
		{	int 
				p=3*i,
				grey = v[p + 2];
			if (grey < 0)grey = 0;
			if (grey >= nshades)	grey = nshades-1;
			if (i*nshades+grey != lg) 
			{	lg = i*nshades+grey;
				g.setColor(col[grey][0]);
			}
			g.fillOval
			(	v[p], v[p + 1],
				2, 2 // SIZE OF THE POINT SOURCE
			);
		}
	}
	/** Find the bounding box of this model */
	void findBB() {
	if (nvert <= 0)
		return;
	float v[] = vert;
	float xmin = v[0], xmax = xmin;
	float ymin = v[1], ymax = ymin;
	float zmin = v[2], zmax = zmin;
	for (int i = nvert * 3; (i -= 3) > 0;) 
	{	float x = v[i];
		if (x < xmin)
		xmin = x;
		if (x > xmax)
		xmax = x;
		float y = v[i + 1];
		if (y < ymin)
		ymin = y;
		if (y > ymax)
		ymax = y;
		float z = v[i + 2];
		if (z < zmin)
		zmin = z;
		if (z > zmax)
		zmax = z;
	}
	this.xmax = xmax;
	this.xmin = xmin;
	this.ymax = ymax;
	this.ymin = ymin;
	this.zmax = zmax;
	this.zmin = zmin;
	}
	//PARTICLE ROUTINES:
	static public void	setVair(double degrees)
	{	double vel=1.0,
			radians=Math.PI*degrees/180.0;//Velocity magnitude
		vdir=degrees;
		Vair[0]=vel*Math.cos(radians);
		Vair[1]=vel*Math.sin(radians);
		Vair[2]=0.0;
	}
	static public double[]	getVair()
	{
		return Vair;
	}
	static public void	setXinj(double xinj, double yinj, double zinj)
	{	Xinj[0]=xinj;
		Xinj[1]=yinj;
		Xinj[2]=zinj;
		for(int i=0;i<DIM;i++)
		{	float d=0.1F;///1.0F/(float)(Dim[i]-2);
			if(Xinj[i]<d)
			{	if(Xinj[i]<0.0)System.out.println("WARNING: Injectction coordinate "+Xinj[i]+" should be positive. Set to 0.0.\n");
				Xinj[i]=d;
			}
			if(Xinj[i]>1.0F-d)
			{	if(Xinj[i]>1.0)System.out.println("WARNING: Injectction coordinate "+Xinj[i]+" should be no greater than 1. Set to 1.\n");
				Xinj[i]=1.0F-d;
			}
		}
	}
	static public void initinj
	(	int n, 
		int red, 
		int green, 
		int blue, 
		float xinj, 
		float yinj, 
		float zinj,
		float vdir
	)
	{	if(n<0) n=0;
		Vair = new double[DIM];
		Vair[0]=1.0;
		Vair[1]=0.0;
		Vair[2]=0.0;
		minj=n;
		ninj=minj;
		Cinj = new int[DIM];
		Cinj[0]=red;
		Cinj[1]=green;
		Cinj[2]=blue;
		for(int i=0;i<3;i++)
		{
			if(Cinj[i]<0)
			{	System.out.println
				("WARNING: Injectction color "+Cinj[i]+" should be positive. Set to 0.\n");
				Cinj[i]=0;
			}
			if(Cinj[i]>255)
			{	System.out.println
				("WARNING: Injectction color "+Cinj[i]+" should be no greater than 255. Set to 255.\n");
				Cinj[i]=255;
			}
		}
		Xinj = new double[DIM];
		setXinj((double)xinj,(double)yinj,(double)zinj);
		Vair = new double[DIM];
		setVair(vdir);
	}
	public int getninj()
	{
		return ninj;
	}
	public void settrapped(int n)
	{
		trapped=n;
	}
	public int gettrapped()
	{
		return trapped;
	}
	public int[] getcolinj()
	{
		return Cinj;
	}
	static public double[] getXinj()
	{ return Xinj;
	}
	static public double getinj(int i)
	{	return Xinj[i];
	}
	static public double getvdir()
	{	return vdir;
	}
	public void setinj()
	{	trapped=0;
		nvert=mvert;
		ninj=minj;
		//SET INJECTION COORDINATES
		for(int i=inj0;i<nvert;i++)
		///for(int i=inj1;i<nvert;i++)
		{	int ip=DIM*i;
			for (int j=0;j<3;j++)
				vert[ip+j]=(float)Math.floor(DX[j]*Xinj[j]);
		}
		int jnj0=DIM*inj0,jnj1=DIM*inj1; 
		for (int i=0;i<DIM;i++)
		{	vert[jnj1+i]=vert[jnj0+i]+vairscale*(float)Vair[i];
		}
	}
///	public int runinj(int nobj, DomainObject Objs[])
	public int runinj(int nobj, int hits[])
	{	int trapped=0,
			i0=nvert-ninj+2,
			ii[] = new int[DIM];
		for(int iobj=0;iobj<nobj;iobj++)
			hits[iobj]=0;
///			Obj[iobj].hits=0;
		float 
			turb=2.0F, // corrected from 2.0 to match the C++
			dt=0.5F;
		boolean outside=false;
		do
		{	outside=false;
			for(int i=i0;i<nvert;i++)
			{	int ip=DIM*i;
				for (int j=0;j<DIM;j++)
				{	int k=(int)((vert[ip+j])/dX[j]+0.5);
					if(k<0) {k=0; outside=true; break;}
					else if(k>=Dim[j]){k=Dim[j]-1;outside=true;break;}
					ii[j]=k;
				}
				///if(outside){trapped++; continue;}
				if(outside)
					if(nvert>delVert(i))
					{	ninj--; i0=i; 
						break;
					}
					else 
					{ 	trapped++;
///						Objs[0].hits++;
						hits[0]++;
						outside=false;
						continue;
					}
				if((int)grid[Domain.I(ii[0],ii[1],ii[2])]!=0)
				{	int iobj=(int)grid[Domain.I(ii[0],ii[1],ii[2])];
					if(iobj<nobj) hits[iobj]++; ///Objs[iobj].hits++;
					trapped++; 
					continue;
				}
				for (int j=0;j<DIM;j++)
				{	float vel=(float)Vair[j]+(float)turb*(-1.0F+2.0F*(float)Math.random());
					vert[ip+j]+=vel*dt;
				}
			}
		}	while(outside);
		return trapped;
	}
}

/** The atomized suface */

/** The representation of a 3D suface */
//+class BoundarySurface
//+{
//+	float vert[];
//+	Atom atoms[];
//+	int tvert[];
//+	int ZsortMap[];
//+	int nvert, maxvert;
//+	
//+	static Hashtable atomTable = new Hashtable();
//+	static Atom defaultAtom;
//+	static 
//+	{
//+		defaultAtom = new Atom(255, 100, 200);
//+	}
//+	boolean transformed;
//+	Matrix3D mat;
//+	
//+	float xmin, xmax, ymin, ymax, zmin, zmax;
//+	
//+	BoundarySurface() 
//+	{
//+		mat = new Matrix3D();
//+		mat.xrot(20);
//+		mat.yrot(30);
//+	}
//+	// Create Atom table
//+	static public void initAtoms (Space space) 
//+	{
//+		atomTable=null;
//+		// Initialize atom table
//+		atomTable = new Hashtable();
//+		for(int icomp=0;icomp<space.vComponents.size();icomp++)
//+		{	Components	comp=(Components)space.vComponents.elementAt(icomp);
//+			atomTable.put
//+			(	comp.name, 
//+				new Atom
//+				(	comp.color.getRed(), 
//+					comp.color.getGreen(), 
//+					comp.color.getBlue()
//+				)
//+			);
//+		}
//+	}
//+	/** Create a boundary surface */
//+	public void init 
//+	(	Space	space,
//+		Domain	domain
//+	) throws Exception
//+	{
//+		double
//+			dx[] = space.getdX(),
//+			x0[] = domain.getX0();
//+		vert=null;
//+		tvert=null;
//+		ZsortMap=null;
//+		nvert=0;
//+		for(int icomp=1;icomp<space.vComponents.size();icomp++)
//+		{	Components	comp=(Components)space.vComponents.elementAt(icomp);
//+			if(comp.boundaryCell!=null)
//+			{	int[] ind=comp.boundaryCell.getInd();
//+				for
//+				(	BoundaryCell curr=comp.boundaryCell;
//+					curr!=null; curr=curr.next
//+				)
//+				{	double
//+						x=x0[0]+curr.ind[0]*dx[0],
//+						y=x0[1]+curr.ind[1]*dx[1],
//+						z=x0[2]+curr.ind[2]*dx[2];
//+					addVert(comp.name, (float) x, (float) y, (float) z);
//+				}
//+			}
//+			else
//+				System.out.println("No boundary for component "+icomp);
//+		}
//+	}  // end BoundarySurface ()
//+	/** Add a vertex to this model */
//+	int addVert(String name, float x, float y, float z) 
//+	{	int i = nvert;
//+		if (i >= maxvert)
//+		if (vert == null) 
//+		{	maxvert = 100;
//+			vert = new float[maxvert * 3];
//+			atoms = new Atom[maxvert];
//+		} else 
//+		{	maxvert *= 2;
//+			float nv[] = new float[maxvert * 3];
//+			System.arraycopy(vert, 0, nv, 0, vert.length);
//+			vert = nv;
//+			Atom na[] = new Atom[maxvert];
//+			System.arraycopy(atoms, 0, na, 0, atoms.length);
//+			atoms = na;
//+		}
//+		Atom a = (Atom) atomTable.get(name);
//+		if (a == null) a = defaultAtom;
//+		atoms[i] = a;
//+		i *= 3;
//+		vert[i] = x;
//+		vert[i + 1] = y;
//+		vert[i + 2] = z;
//+		return nvert++;
//+	}
//+	/** Transform all the points in this model */
//+	void transform() 
//+	{	if (transformed || nvert <= 0) return;
//+		if (tvert == null || tvert.length < nvert * 3)
//+		tvert = new int[nvert * 3];
//+		mat.transform(vert, tvert, nvert);
//+		transformed = true;
//+	}
//+	/** Paint this model to a graphics context.  It uses the matrix associated
//+	with this model to map from model space to screen space.
//+	The next version of the browser should have double buffering,
//+	which will make this *much* nicer */
//+	void paint(Graphics g) 
//+	{
//+		if (vert == null || nvert <= 0) return;
//+		transform();
//+		int v[] = tvert;
//+		int zs[] = ZsortMap;
//+		if (zs == null) 
//+		{
//+			ZsortMap = zs = new int[nvert];
//+			for (int i = nvert; --i >= 0;)
//+			zs[i] = i * 3;
//+		}
//+		/*
//+		* I use a bubble sort since from one iteration to the next, the sort
//+		* order is pretty stable, so I just use what I had last time as a
//+		* "guess" of the sorted order.  With luck, this reduces O(N log N)
//+		* to O(N)
//+		*/
//+		for (int i = nvert - 1; --i >= 0;) 
//+		{	boolean flipped = false;
//+			for (int j = 0; j <= i; j++) 
//+			{	int a = zs[j];
//+				int b = zs[j + 1];
//+				if (v[a + 2] > v[b + 2]) 
//+				{
//+					zs[j + 1] = a;
//+					zs[j] = b;
//+					flipped = true;
//+				}
//+			}
//+			if (!flipped)
//+			break;
//+		}
//+		int lg = 0;
//+		int lim = nvert;
//+		Atom ls[] = atoms;
//+		if (lim <= 0 || nvert <= 0)return;
//+		for (int i = 0; i < lim; i++) 
//+		{	int j = zs[i];
//+			int grey = v[j + 2];
//+			if (grey < 0)
//+			grey = 0;
//+			if (grey > 15)
//+			grey = 15;
//+			// g.drawString(names[i], v[j], v[j+1]);
//+			atoms[j/3].paint(g, v[j], v[j + 1], grey);
//+			// g.drawImage(iBall, v[j] - (iBall.width >> 1), v[j + 1] -
//+			// (iBall.height >> 1));
//+		}
//+	}
//+	/** Find the bounding box of this model */
//+	void findBB() 
//+	{	if (nvert <= 0)
//+		return;
//+		float v[] = vert;
//+		float xmin = v[0], xmax = xmin;
//+		float ymin = v[1], ymax = ymin;
//+		float zmin = v[2], zmax = zmin;
//+		for (int i = nvert * 3; (i -= 3) > 0;) 
//+		{	float x = v[i];
//+			if (x < xmin)
//+			xmin = x;
//+			if (x > xmax)
//+			xmax = x;
//+			float y = v[i + 1];
//+			if (y < ymin)
//+			ymin = y;
//+			if (y > ymax)
//+			ymax = y;
//+			float z = v[i + 2];
//+			if (z < zmin)
//+			zmin = z;
//+			if (z > zmax)
//+			zmax = z;
//+		}
//+		this.xmax = xmax;
//+		this.xmin = xmin;
//+		this.ymax = ymax;
//+		this.ymin = ymin;
//+		this.zmax = zmax;
//+		this.zmin = zmin;
//+	}
//+}


/** An applet to put a model into a page */
public class DomainView
//extends JPanel
extends Applet
implements Runnable, MouseListener, MouseMotionListener 
{
//+	static final int MODE_WIREFRAME = 0;
//+	static final int MODE_ATOMS = 1;
//+	int mode = MODE_WIREFRAME;

	static Space	space;
	static Domain	domain=null;
//-	BoundarySurface md;
	boolean painted = true;
	boolean doRun = false;       // True if thread currently paused by user
	float xfac;
	int prevx, prevy;
	float xtheta, ytheta;
	float scalefudge = 1;
	Matrix3D amat = new Matrix3D(), tmat = new Matrix3D();
	String mdname = null;
	String message = null;
	Image backBuffer;
	Graphics backGC = null;
	Dimension backSize;
	WireFrame wireframe;
///	JComboBox comboBox3DRenderModes = new JComboBox();

	float deltaX = 0.0f, deltaY = 0.0f, scale = 1.0f;
	boolean left = false, center = false, right = false;
    
	JPanel drawingSurface = new JPanel();
///    JButton buttonRotX = new JButton("Rotate X"), buttonRotY = new JButton("Rotate Y");
	JScrollBar scrollBarX,scrollBarY,scrollBarZ,scrollBarA;
	JButton buttonRun = new JButton("Fire!");
	JPanel dataPanel;
///	JPanel dataTable;
///	JLabel dataLabel[];
	JLabel 
	 	objectColorLabel,
		predictedHitsTitle,predictedHitsLabel,
		actualHitsTitle,actualHitsLabel;
//?	JComboBox<String> selectObjectComboBox;
	JComboBox selectObjectComboBox;
	//Shooting area plane indicator

	int nobjs=0, //max number of objs
		hits[]=null,
		nobjx, nobjy, nobjz, nobja,
		Hits[]=null,
		iobj=0; //obj of the object to count hits for
	String objname[];
	Color objcolor[];
///	DomainObject Objs[];
	
	private synchronized void newBackBuffer() 
	{
		backBuffer = createImage(getSize().width, getSize().height);
		if (backGC != null) 
		{	backGC.dispose();
		}
		backGC = backBuffer.getGraphics();
		backSize = getSize();
	}
	public class SelectObject implements ActionListener
	{	public void actionPerformed(ActionEvent evt)
		{	iobj = selectObjectComboBox.getSelectedIndex();
			if(iobj>0)
			{	//LOAD SAMPLE FILE:
				///ObjData.Load(iobj, getDocumentBase());
///				boolean outside=false;
				int	DIM=space.getDim(), 
					IX[]=new int[DIM],
					X0[]=ObjData.getX0(),
					Dim[]=ObjData.getDim();
				double 
					dY[] = space.getdX(),
					DY[] = domain.getDX(),
					y[] = wireframe.getXinj();
				for(int i=0;i<DIM;i++)
				{	//Cell-index of the point
					IX[i]=(int)Math.floor(DY[i]*y[i]/dY[i]);
///					if(IX[i]<X0[i]||IX[i]>=X0[i]+Dim[i])
///						outside=true;//data not available
				}
				// Load data
				if(ObjData.Load(iobj,IX,getDocumentBase())>=0)
				{	wireframe.setFrame();
					objectColorLabel.setText
					("Object: "+iobj+", Data Set: "+ObjData.getDataIndex());
				}
				else
				{	wireframe.hideFrame();
					objectColorLabel.setText("Object: "+iobj+", Data Set: N/A");
				}
				objectColorLabel.setOpaque(true);
				repaint();
			}// end LOAD SAMPLE FILE
			else 
			{	Hits=null;
				objectColorLabel.setOpaque(false);
			}
			objectColorLabel.setBackground(objcolor[iobj]);
		}
	}
	public void init() 
	{
///		comboBox3DRenderModes.addItem("Wireframe View");
///		comboBox3DRenderModes.addItem("Atom View");
///		comboBox3DRenderModes.addActionListener(new ChangeRenderMode());
///		mode = MODE_WIREFRAME;
		int ninj=0,//number of injected particles
			rinj=255, //particles default colors
			ginj=0,
			binj=0;
		float 
			xrot0=0.0F,
			yrot0=0.0F,
			zrot0=0.0F,
			xinj=0.0F,
			yinj=0.0F,
			zinj=0.0F,
			vdir=0.0F;
		String str=null;
		mdname = getParameter("model");
		if (mdname == null)
			mdname = "project.config";
		try
		{	str = getParameter("xrot");
			if (str != null) 
				xrot0 = (Float.valueOf(str)).floatValue() ;
			str = getParameter("yrot");
			if (str != null) 
				yrot0 = (Float.valueOf(str)).floatValue() ;
			str = getParameter("zrot");
			if (str != null) 
				zrot0 = (Float.valueOf(str)).floatValue() ;
			str = getParameter("ninj");
			if (str != null) 
				ninj = (Integer.valueOf(str)).intValue() ;
			str = getParameter("xinj");
			if (str != null) 
				xinj = (Float.valueOf(str)).floatValue() ;
			if(xinj<0.0)xinj=0.0F; if(xinj>1.0)xinj=1.0F;
			str = getParameter("yinj");
			if (str != null) 
				yinj = (Float.valueOf(str)).floatValue() ;
			if(yinj<0.0)yinj=0.0F; if(yinj>1.0)yinj=1.0F;
			str = getParameter("zinj");
			if (str != null) 
				zinj = (Float.valueOf(str)).floatValue() ;
			if(zinj<0.0)zinj=0.0F; if(zinj>1.0)zinj=1.0F;
			str = getParameter("rinj");
			if (str != null) 
				rinj = (Integer.valueOf(str)).intValue() ;
			if(rinj<0)rinj=0; if(rinj>255)rinj=255;
			str = getParameter("ginj");
			if (str != null) 
				ginj = (Integer.valueOf(str)).intValue() ;
			if(ginj<0)ginj=0; if(ginj>255)ginj=255;
			str = getParameter("binj");
			if (str != null) 
				binj = (Integer.valueOf(str)).intValue() ;
			if(binj<0)binj=0; if(binj>255)binj=255;
			str = getParameter("vdir");
			if (str != null) 
				vdir = (Float.valueOf(str)).floatValue() ;
		}
		catch (NumberFormatException e) 
		{///	System.out.println("Failed to read from "+mdname+": "+e.toString());
			System.out.println("Parameter error: "+e.toString());
		}
		String	projectDirectory = "./";
		space= new Space("Demo",1.0,1.0,1.0);
		URL url = this.getCodeBase();
		String url_str = url.toString();
///		System.out.println("URL: "+url_str);
		space.setPath(projectDirectory);
		space.vComponents.clear();
		Components empty = new Components();
		empty.name = "This is an empty component.  It's here only to make things compatible w/the master.";
		empty.color = new Color(255,255,255);
		space.vComponents.add(empty);
		try
		{	InputStream is = new URL(getDocumentBase(),mdname).openStream();
			StreamTokenizer in = new StreamTokenizer
			(new BufferedReader(new InputStreamReader(is)));
			//Read in the components
			in.nextToken(); //NumComponents
			in.nextToken();
			int numComponents = (int)in.nval;
			nobjs=numComponents+1;
///			Objs = new DomainObject[nobjs];
			hits=new int[nobjs];
			objname=new String[nobjs];
			objcolor=new Color[nobjs];
			objname[0]="SELECT OBJECT";
			objcolor[0]=new Color(255,255,255);
			for ( int i = 0; i < numComponents; i++ )
			{	Components newDCI = new Components();
				in.nextToken();
				newDCI.name = in.sval;
///				Objs[i].name=new String(newDCI.name);
				objname[i+1]=new String(newDCI.name);
				hits[i+1]=0;
				in.nextToken();
				int r = (int)in.nval;
				in.nextToken();
				int g = (int)in.nval;
				in.nextToken();
				int b = (int)in.nval;
				newDCI.color = new Color(r, g, b);
///				Objs[i].color=new Color(newDCI.color);
				objcolor[i+1]=new Color(r,g,b);
				space.vComponents.add(newDCI);
			}
			///if (space.vComponents.size() > 1)
			//Read in the domain names and load them as "ghosts"
			in.nextToken(); //NumDomains
			in.nextToken();
			int numDomains = (int)in.nval;///DDD: this should be removed
			                             /// we don't need the domain number
			//Empty all the current domains/connections
			space.RemoveDomains();
			Domain dom = new Domain("dummy", 0, 0, 0, space);
			for (int index = 0; index < numDomains; index++)
			{	in.nextToken();
				/*
				File f = new File
				(projectDirectory + File.separator + in.sval + ".dom");
				*/
				//if (f.exists())
				{	dom.name = in.sval;
					dom.LoadGhost(url_str);
					if (space.first_domain == null)
						space.first_domain = dom;
					else
						space.first_domain.Append(dom);
				}
				/*
				else
				JOptionPane.showMessageDialog
				(	null,
					new JLabel("Domain " + in.sval + " doesn't exist.  Did you delete it?"), 
					"Warning!",
					JOptionPane.WARNING_MESSAGE
				);
				*/
			}
			setLayout(new BorderLayout());
///			drawingSurface.setVisible(false);
			add(drawingSurface);
			JPanel controlPanel = new JPanel();
///			buttonRotX.addActionListener(new Rotate());
///			buttonRotY.addActionListener(new Rotate());
///			controlPanel.add(comboBox3DRenderModes);
///			controlPanel.add(buttonRotX); 
///			controlPanel.add(buttonRotY);

			scrollBarX = new JScrollBar
			(	JScrollBar.VERTICAL, 1, 1, 1, 
				dom.Dim[0]-1
			);
			scrollBarX.addAdjustmentListener(new ChangeXinj());
			scrollBarX.setToolTipText("Injection position X");
			controlPanel.add(scrollBarX); 

			scrollBarY = new JScrollBar
			(	JScrollBar.VERTICAL, 1, 1, 1, 
				dom.Dim[1]-1
			);
			scrollBarY.addAdjustmentListener(new ChangeYinj());
			scrollBarY.setToolTipText("Injection position Y");
			controlPanel.add(scrollBarY); 

			scrollBarZ = new JScrollBar
			(	JScrollBar.VERTICAL, 1, 1, 1, 
				dom.Dim[2]-1
			);
			scrollBarZ.addAdjustmentListener(new ChangeZinj());
			scrollBarZ.setToolTipText("Injection position Z");
			controlPanel.add(scrollBarZ); 

			scrollBarA = new JScrollBar
			(	JScrollBar.VERTICAL, 1, 1, 1, 
				360
			);
			scrollBarA.addAdjustmentListener(new ChangeVair());
			scrollBarA.setToolTipText("Wind direction angle");
			controlPanel.add(scrollBarA); 

			buttonRun.addActionListener(new Run());
			controlPanel.add(buttonRun); 

			add(controlPanel);

			dataPanel = new JPanel(new GridLayout(3,2));
			dataPanel.setAlignmentX(JPanel.LEFT_ALIGNMENT);
///			dataPanel.setOpaque(false);
///			dataPanel.setBackground(new Color(255,255,255));
			// 1,1:
//?			selectObjectComboBox= new JComboBox<String>(); 
			selectObjectComboBox= new JComboBox(); 
			for(int i=0;i<nobjs;i++)
			{	String color[]=new String[3];
				color[0]=Integer.toHexString(objcolor[i].getRed());
				color[1]=Integer.toHexString(objcolor[i].getGreen());
				color[2]=Integer.toHexString(objcolor[i].getBlue());
				for(int j=0;j<3;j++) if(color[j].length()<2) color[j]="0"+color[j];
				selectObjectComboBox.addItem
				(	"<html><font color=\"#"
					+color[0]
					+color[1]
					+color[2]
					+"\"><b>***</b></font> "
					+objname[i]
					+"</html>"
				); 
			}
///				selectObjectComboBox.addItem(Objs[i].name); 
			selectObjectComboBox.addActionListener(new SelectObject());
			selectObjectComboBox.setAlignmentX(JComboBox.LEFT_ALIGNMENT); 
			dataPanel.add(selectObjectComboBox);  
			// 1,2:
			objectColorLabel = new JLabel();
///			objectColorLabel.setOpaque(true);
			dataPanel.add(objectColorLabel);
			//2,1:
			predictedHitsTitle = new JLabel("Predicted:");
			dataPanel.add(predictedHitsTitle);
			//2,2:
			predictedHitsLabel = new JLabel("0");
			dataPanel.add(predictedHitsLabel);
			//2,1:
			actualHitsTitle = new JLabel("Current:");
			dataPanel.add(actualHitsTitle);
			//3,2:
			actualHitsLabel = new JLabel("0");
			dataPanel.add(actualHitsLabel);
			
			add(dataPanel, BorderLayout.SOUTH);

///			actualHitsLabel.repaint();
///			dataPanel.setVisible(true);
///			hits = new int[nobjs];
		}
		catch(Exception e)
		{
			System.out.println("exception in init():"+ e);
		}
		try 
		{	str = getParameter("scale");
			if(str!=null)scalefudge = Float.valueOf(str).floatValue();
		}	catch(Exception e) 
		{
			System.out.println("exception in init() valueOf:"+ e);
		}
		amat.xrot(xrot0);
		amat.yrot(yrot0);
		amat.zrot(zrot0);
		resize(getSize().width <= 300 ? 480 : getSize().width,
		getSize().height <= 300 ? 480 : getSize().height);
		newBackBuffer();
		addMouseListener(this);
		addMouseMotionListener(this);
		domain= space.first_domain;
		domain.Load(url_str);
		domain.setBoundaryCells
		(
			space.vComponents
		);
		for 
		(	int plane_orientation = 0;
				plane_orientation < space.getDim();
			plane_orientation++
		)
		{	int separation=1;///DDD should be user-defined via gui.
			domain.setBoundaryContours
			(
				plane_orientation, ///DDD: should be user-defined
				separation, ///DDD: should be user-defined
				space.vComponents
			);
		}
//+		BoundarySurface.initAtoms(space);
		WireFrame.initinj(ninj,rinj,ginj,binj,xinj,yinj,zinj,vdir);
		wireframe = new WireFrame (space, domain, space.vComponents);
		scrollBarA.setValue((int)(WireFrame.getvdir()));
		scrollBarX.setValue((int)(WireFrame.getinj(0)*(double)domain.Dim[0]));
		scrollBarY.setValue((int)(WireFrame.getinj(1)*(double)domain.Dim[1]));
		scrollBarZ.setValue((int)(WireFrame.getinj(2)*(double)domain.Dim[2]));
		validate();

/// MOVED FROM run():
		try 
		{	Thread.currentThread().setPriority(Thread.MIN_PRIORITY);
		//if (domain==null)break;
///			WireFrame m = new WireFrame (space, domain, space.vComponents);
///			wireframe = m;
			WireFrame m = wireframe;
			m.findBB();
			m.compress(); ///QS
			float xw = m.xmax - m.xmin;
			float yw = m.ymax - m.ymin;
			float zw = m.zmax - m.zmin;
			if (yw > xw)
			xw = yw;
			if (zw > xw)
			xw = zw;
			float f1 = getSize().width / xw;
			float f2 = getSize().height / xw;
			xfac = 0.7f * (f1 < f2 ? f1 : f2) * scalefudge;

		} 
		catch(Exception e) 
		{
			wireframe = null;
			message = e.toString();
		}
///		run();
//+		break;
//+		}
/// MOVED FROM run

		repaint();
	}
	static public void init(Space space, Domain D)
	{
		domain=D;
//+		BoundarySurface.initAtoms(space);
	}
	public void destroy() 
	{
		removeMouseListener(this);
		removeMouseMotionListener(this);
	}
	public void run() 
	{	int
			mp=wireframe.getninj()-2,//number of free particles (excludes the injection point)
			np=mp,
			trapped=wireframe.gettrapped();//number of trapped particles
			do //RUN PARTICLES
			{	if(doRun)	
				{	trapped=wireframe.runinj(nobjs,hits);///,Objs);
					np=wireframe.getninj()-2; 
					repaint();
					if(iobj<nobjs)
//					actualHitsLabel.setText(Integer.toString(hits[iobj]));
					{	int error=100,
							probability=(int)(100.0F*(float)hits[iobj]/(float)mp);
						if(hits[iobj]>0) 
							error=(int)Math.floor(200.0/Math.sqrt((float)hits[iobj]));
						String errorlabel=error<1?"< 1":Integer.toString(error);
						if(probability<=0) errorlabel="N/A";
						actualHitsLabel.setText
						(	Integer.toString(probability)+"%"
							+" Error: "+errorlabel+"%"
						);
					}
				}
			}	while(doRun&&trapped<np);
		doRun=false;
	}
	public void start() 
	{
//-	   if (md == null && message == null)
	   if (message == null)
			new Thread(this).start();
	}
	public void stop() 
	{
	}
	/* event handling */
	public void mouseClicked(MouseEvent e) 
	{
	}
	public void mousePressed(MouseEvent e) 
	{
		if (e.getButton() == e.BUTTON1)
		left = true;
		else if (e.getButton() == e.BUTTON2)
		center = true;
		else if (e.getButton() == e.BUTTON3)
		right = true;
		prevx = e.getX();
		prevy = e.getY();
		e.consume();
	}
	public void mouseReleased(MouseEvent e) 
	{
		if (e.getButton() == e.BUTTON1)
		left = false;
		else if (e.getButton() == e.BUTTON2)
		center = false;
		else if (e.getButton() == e.BUTTON3)
		right = false;
	}
	public void mouseEntered(MouseEvent e) 
	{
	}
	public void mouseExited(MouseEvent e) 
	{	
	}
	public void mouseDragged(MouseEvent e) 
	{	int x = e.getX();
		int y = e.getY();

		if (center || (left && right))
		{
			deltaX += -(prevx - x);
			deltaY += -(prevy - y);
		}
		else if (left)
		{
			tmat.unit();
			float xtheta = (prevy - y) * (360.0f / getSize().width);
			float ytheta = (x - prevx) * (360.0f / getSize().height);
			tmat.xrot(xtheta);
			tmat.yrot(ytheta);
			amat.mult(tmat);
		}
		else if (right)
		{
			scale += (prevy - y)/180.0f;
			if (scale < 0.1f)
			    scale = 0.1f;
		}
		if (painted) 
		{	painted = false;
			repaint();
		}
		prevx = x;
		prevy = y;
		e.consume();
	}
	public void mouseMoved(MouseEvent e) 
	{
	}
	public void update(Graphics g) 
	{	if (backBuffer == null)
		g.clearRect(0, 0, getSize().width, getSize().height);
		paint(g);
	}
	public void paint(Graphics g) 
	{
	    //g.clearRect(0,0,getSize().width, getSize().height); //This is ENTIRELY too slow

//+		switch(mode)
//+		{
//+		case MODE_ATOMS:
//+		if (md != null) 
//+		{	md.mat.unit();
//+			md.mat.translate(-(md.xmin + md.xmax) / 2,
//+			-(md.ymin + md.ymax) / 2,
//+			-(md.zmin + md.zmax) / 2);
//+			md.mat.mult(amat);
//+			md.mat.scale(xfac, -xfac, 16 * xfac / getSize().width);
//+			md.mat.scale(scale,scale,scale);
//+			md.mat.translate(getSize().width / 2, getSize().height / 2, 8);
//+			md.mat.translate(deltaX, deltaY, 0.0f);
//+			md.transformed = false;
//+			if (backBuffer != null) 
//+			{	if (!backSize.equals(getSize()))
//+				newBackBuffer();
//+				backGC.setColor(getBackground());
//+				backGC.fillRect(0,0,getSize().width,getSize().height);
//+				md.paint(backGC);
//+				g.drawImage(backBuffer, 0, 0, this);
//+			}
//+			else
//+				md.paint(g);
//+			setPainted();
//+		} else if (message != null) 
//+		{	g.drawString("Error in model:", 3, 20);
//+			g.drawString(message, 10, 40);
//+		}
//+		break;
//+		case MODE_WIREFRAME: // WireFrame
		if (wireframe != null) 
		{
			//Paint scene
			wireframe.mat.unit();
			wireframe.mat.translate(-(wireframe.xmin + wireframe.xmax) / 2,
					-(wireframe.ymin + wireframe.ymax) / 2,
					-(wireframe.zmin + wireframe.zmax) / 2);
			wireframe.mat.mult(amat);
			wireframe.mat.scale(xfac, -xfac, 16 * xfac / getSize().width);
			wireframe.mat.scale(scale,scale,scale);
			wireframe.mat.translate(getSize().width / 2, getSize().height / 2, 8);
			wireframe.mat.translate(deltaX, deltaY, 0.0f);
			wireframe.transformed = false;
			if (backBuffer != null)
			{
				if (!backSize.equals(getSize()))
					newBackBuffer();
				backGC.setColor(getBackground());
				backGC.fillRect(0,0,getSize().width, getSize().height);
				//Clear the screen (very slow)
				wireframe.paint(backGC);
				g.drawImage(backBuffer, 0, 0, this);
			}
			else
				wireframe.paint(g);
			setPainted();
		}
		else if (message != null) 
		{	g.drawString("Error in model:", 3, 20);
			g.drawString(message, 10, 40);
		}
//+		break;
//+		}
///		comboBox3DRenderModes.repaint(); 
///		buttonRotX.repaint();
///		buttonRotY.repaint();
		buttonRun.repaint();
		scrollBarA.repaint();
		scrollBarX.repaint();
		scrollBarY.repaint();
		scrollBarZ.repaint();

		objectColorLabel.repaint();
		selectObjectComboBox.repaint();
		actualHitsTitle.repaint();
		actualHitsLabel.repaint();
		predictedHitsTitle.repaint();
		predictedHitsLabel.repaint();
	}
	private synchronized void setPainted() 
	{	painted = true;
		notifyAll();
	}
	private synchronized void waitPainted()
	{	while (!painted)
		{	try
			{
				wait();
			}
			catch (InterruptedException e) 
			{
			}
		}
		painted = false;
	}
	public String getAppletInfo() 
	{
		return "Title: DomainView \nAuthor: A.Smirnov, B.Sowers, H.Zhang \nDomain 3D Viewer.";
	}
	public String[][] getParameterInfo() 
	{	String[][] info = 
		{	{"model", "path string", "The path to the model to be displayed in .xyz format (see http://chem.leeds.ac.uk/Project/MIME.html).  Default is model.obj."
			},
			{"scale", "float", "Scale factor.  Default is 1 (i.e. no scale)."
			}
		};
		return info;
	}
	public void rotx(float angle)
	{
		tmat.unit();
		tmat.xrot(angle);
		amat.mult(tmat);
		repaint();
	}
	public void roty(float angle)
	{
		tmat.unit();
		tmat.yrot(angle);
		amat.mult(tmat);
		repaint();
	}
	class Run implements ActionListener
	{	public void actionPerformed(ActionEvent evt)
		{	if (evt.getSource() == buttonRun) 
			{	wireframe.setinj();
				// Check if the data set needs to be reloaded
				boolean outside=false;
				int	DIM=space.getDim(), 
					IX[]=new int[DIM],
					X0[]=ObjData.getX0(),
					Dim[]=ObjData.getDim();
				double 
					dY[] = space.getdX(),
					DY[] = domain.getDX(),
					y[] = wireframe.getXinj();
				for(int i=0;i<DIM;i++)
				{	//Cell-index of the point
					IX[i]=(int)Math.floor(DY[i]*y[i]/dY[i]);
					if(IX[i]<X0[i]||IX[i]>=X0[i]+Dim[i])
						outside=true;//data not available
				}
				// Reload data if needed
				int idata=ObjData.getDataIndex();
				if(outside) 
				{	idata=ObjData.Load(iobj,IX,getDocumentBase());
				}
				if(idata>=0)
				{	wireframe.setFrame();
					objectColorLabel.setText
					("Object: "+iobj+", Data Set: "+ObjData.getDataIndex());
				}
				else
				{	wireframe.hideFrame();
					objectColorLabel.setText("Object: "+iobj+", Data Set: N/A");
				}
				objectColorLabel.setOpaque(true);
				Statistics hitStatistics=new Statistics();	
				ObjData.getHitStatistics
				(	space.getDim(),space.getdX(),domain.getDX(),
					wireframe.getXinj(),wireframe.getvdir(),
					hitStatistics
				);
//				float hitProbability=ObjData.getHitProbability
//				(	space.getDim(),space.getdX(),domain.getDX(),wireframe.getXinj(),wireframe.getvdir()
//				);
				if(hitStatistics.getProbability()>=0.0)
				{	float p=hitStatistics.getProbability();
					int probability=(int)Math.round(100.0F*p),
						error=100;
					String errorlabel="N/A";
					if(probability>0) 
					{	error=(int)Math.round(100.0F*hitStatistics.getError());
						errorlabel=error<1?"< 1":Integer.toString(error);
					}
					predictedHitsLabel.setText
					(	Integer.toString(probability)+"%"
						+" Error: "+errorlabel+"%"
					);
				}
				else
					predictedHitsLabel.setText("N/A");
				doRun=true;	
				start();
//				repaint();
			}
		}
	}
	class ChangeVair implements AdjustmentListener
	{	public void adjustmentValueChanged(AdjustmentEvent e)
		{	if (WireFrame.getVair()!=null)
			{	doRun=false;
				WireFrame.setVair(e.getValue());
				wireframe.setinj();
				repaint();
			}
		}
	}
	class ChangeXinj implements AdjustmentListener
	{	public void adjustmentValueChanged(AdjustmentEvent e)
		{	if (WireFrame.getXinj()!=null)
			{	doRun=false;
				WireFrame.setXinj
				(	(double)e.getValue()/(double)(domain.Dim[0]-2),
					WireFrame.getinj(1),
					WireFrame.getinj(2)
				);
				wireframe.setinj();
				repaint();
			}
		}
	}
	class ChangeYinj implements AdjustmentListener
	{	public void adjustmentValueChanged(AdjustmentEvent e)
		{	if (WireFrame.getXinj()!=null)
			{	doRun=false;
				WireFrame.setXinj
				(	WireFrame.getinj(0),
					(double)e.getValue()/(double)(domain.Dim[1]-2),
					WireFrame.getinj(2)
				);
				wireframe.setinj();
				repaint();
			}
		}
	}
	class ChangeZinj implements AdjustmentListener
	{	public void adjustmentValueChanged(AdjustmentEvent e)
		{	if (WireFrame.getXinj()!=null)
			{	doRun=false;
				WireFrame.setXinj
				(	WireFrame.getinj(0),
					WireFrame.getinj(1),
					(double)e.getValue()/(double)(domain.Dim[2]-2)
				);
				wireframe.setinj();
				repaint();
			}
		}
	}

///    public class Rotate implements ActionListener
///    {
///	public void actionPerformed(ActionEvent evt)
///	{
///	    if (evt.getSource() == buttonRotX)
///		rotx(10);
///	    if (evt.getSource() == buttonRotY)
///		roty(10);
///	}
///    }
}   // end class DomainView


///	public class ChangeRenderMode implements ActionListener
///	{
///		public void actionPerformed(ActionEvent evt)
///		{
///			mode = comboBox3DRenderModes.getSelectedIndex();
///			run();
///		}
///	}
//+class Atom 
//+{
//+///	private static JPanel panel;
//+	private static Applet applet;
//+	private static byte[] data;
//+	private final static int R = 40;
//+	private final static int hx = 15;
//+	private final static int hy = 15;
//+	private final static int bgGrey = 192;
//+	private final static int nBalls = 16;
//+	private static int maxr;
//+	
//+	private int Rl;
//+	private int Gl;
//+	private int Bl;
//+	private Image balls[];
//+	
//+	static 
//+	{	data = new byte[R * 2 * R * 2];
//+		int mr = 0;
//+		for (int Y = 2 * R; --Y >= 0;) 
//+		{	int x0 = (int) (Math.sqrt(R * R - (Y - R) * (Y - R)) + 0.5);
//+			int p = Y * (R * 2) + R - x0;
//+			for (int X = -x0; X < x0; X++) 
//+			{	int x = X + hx;
//+				int y = Y - R + hy;
//+				int r = (int) (Math.sqrt(x * x + y * y) + 0.5);
//+				if (r > mr)
//+				mr = r;
//+				data[p++] = r <= 0 ? 1 : (byte) r;
//+			}
//+		}
//+		maxr = mr;
//+	}
//+///	static void setPanel(JPanel p) 
//+///	{
//+///		panel = p;
//+///	}
//+	static void setApplet(Applet app) 
//+	{
//+		applet = app;
//+	}
//+	Atom(int Rl, int Gl, int Bl) 
//+	{
//+		this.Rl = Rl;
//+		this.Gl = Gl;
//+		this.Bl = Bl;
//+	}
//+	private final int blend(int fg, int bg, float fgfactor) 
//+	{
//+		return (int) (bg + (fg - bg) * fgfactor);
//+	}
//+	private void Setup() 
//+	{
//+		balls = new Image[nBalls];
//+		byte red[] = new byte[256];
//+		red[0] = (byte) bgGrey;
//+		byte green[] = new byte[256];
//+		green[0] = (byte) bgGrey;
//+		byte blue[] = new byte[256];
//+		blue[0] = (byte) bgGrey;
//+		for (int r = 0; r < nBalls; r++) 
//+		{	float b = (float) (r+1) / nBalls;
//+			for (int i = maxr; i >= 1; --i) 
//+			{	float d = (float) i / maxr;
//+				red[i] = (byte) blend(blend(Rl, 255, d), bgGrey, b);
//+				green[i] = (byte) blend(blend(Gl, 255, d), bgGrey, b);
//+				blue[i] = (byte) blend(blend(Bl, 255, d), bgGrey, b);
//+			}
//+			IndexColorModel model = new IndexColorModel
//+			(	8, maxr + 1,
//+				red, green, blue, 0
//+			);
//+			balls[r] = applet.createImage
//+			(	new MemoryImageSource
//+				(R*2, R*2, model, data, 0, R*2)
//+			);
//+		}
//+	}
//+	void paint(Graphics gc, int x, int y, int r) 
//+	{
//+		Image ba[] = balls;
//+		if (ba == null) 
//+		{
//+			Setup();
//+			ba = balls;
//+		}
//+		Image i = ba[r];
//+		int size = 10 + r;
//+		gc.drawImage(i, x - (size >> 1), y - (size >> 1), size, size, applet);
//+	}
//+}
