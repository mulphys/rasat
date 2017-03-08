import java.awt.*;
import java.util.*;

public class	BoundaryContour
{
///	static int	DIM, Dim[];
///	static float[]	x0,dx;
	int	ipos; //plane position
	int[]	ind; // plane-indexes of the vertexes: compressed as (ix<<16)|iy
	public BoundaryContour	next;
///	public static void	init
///	(	int	DIM0,
///		int	dim[],
///		double	X0[],
///		double	DX[]
///	)
///	{	DIM=DIM0;
///		Dim = new int[DIM];
///		x0 = new float[DIM];
///		dx = new float[DIM];
///		for (int i=0;i<DIM;i++)
///		{	Dim[i]=dim[i];
///			dx[i]=(float)DX[i]/(float)Dim[i];
///			x0[i]=(float)X0[i]+0.5F*(float)dx[i];
///		}
///	}
	public BoundaryContour(int ipos, int n)
	{
		this.ipos=ipos;
		this.ind = new int[n];
		next=null;
	}
	public int getPlaneInd()
	{
		return ipos;
	}
	public int getVertexIndX(int i)
	{
		return ((ind[i] >> 16) & 0xFFFF);
	}
	public int getVertexIndY(int i)
	{
		return (ind[i] & 0xFFFF);
	}
	public void getVertexIndXY(int i, int[] ind)
	{	int	j=ind[i];
		ind[0] = ((j >> 16) & 0xFFFF);
		ind[1] = (j & 0xFFFF);
	}
	public void putVertex(int i, int val)
	{	try
		{	ind[i]=val;
		}	catch(Exception e)
		{
		    //System.out.println(e.toString());
		}
	}
///	public void setNext(BoundaryContour next) //let's use 'next' directly!
///	{	this.next=next;
///	}
	public int	Length()
	{
		return ind.length;
	}
}

