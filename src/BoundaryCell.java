import java.awt.*;

public class	BoundaryCell
{
	public static int	DIM, Dim[];
	public static float[]	x0,dx;
	public int	ind[]; // grid-indexes of the cell 
	public BoundaryCell	next;
	static public void init
	(	int	DIM0,
		int	dim[],
		double	X0[],
		double	DX[]
	)
	{	DIM=DIM0;
		Dim = new int[DIM];
		x0 = new float[DIM];
		dx = new float[DIM];
		//System.out.println("BoundaryCell: Dim="+dim[0]+","+dim[1]+","+dim[2]);///DDD
		//System.out.println("BoundaryCell:DX="+DX[0]+","+DX[1]+","+DX[2]);///DDD
		for (int i=0;i<DIM;i++)
		{	Dim[i]=dim[i];
			dx[i]=(float)DX[i]/(float)Dim[i];
			x0[i]=(float)X0[i]+0.5F*(float)dx[i];
		}
	}
	public	BoundaryCell()
	{
		ind = new int[DIM];
		next=null;
	}
	public	BoundaryCell
	(
		int	i0,
		int	i1,
		int	i2
	)
	{	this();
		ind[0]=i0;
		ind[1]=i1;
		ind[2]=i2;
	}
	public void setNext(BoundaryCell next)
	{
		this.next=next;
	}
	public float	X(int i)
	{
		return x0[i]+dx[i]*(float)ind[i];
	}
	public int[]	getInd()
	{
		return ind;
	}
}
///class	BoundaryFace
///{
///	static final int	DIM=3;
///	float[]	
///		x, // coordinates of face center
///		v; // unit normal vector to the face
///	BoundaryFace	next;
///	public	BoundaryFace()
///	{
///		x = new float[DIM];
///		v = new float[DIM];
///		next=null;
///	}
///	public	BoundaryFace
///	(
///		float	x0,
///		float	x1,
///		float	x2,
///		float	v0,
///		float	v1,
///		float	v2
///	)
///	{	this();
///		x[0]=x0;
///		x[1]=x1;
///		x[2]=x2;
///		v[0]=v0;
///		v[1]=v1;
///		v[2]=v2;
///	}
///	public void setNext(BoundaryFace next)
///	{
///		this.next=next;
///	}
///}
///class	DomainComponent 
///{
///	static final int	DIM=3;
///	int
///		type, // 0: empty, 1... different component types
///		m0,m01; // used in I(i,j,k)
///	Variable[] vars;
///	BoundaryFace	boundarySurface;
///	DomainComponent	next,prev;
///	public int	I
///	(	int	i0,
///		int	i1,
///		int	i2
///	)
///	{	return i0+i1*m0+i2*m01;
///	}
///	public DomainComponent(int type)
///	{	this.type=type;
///		boundaryCell=null;
///		boundarySurface=null;
///		vars=null;
///	}
///	public void deleteBoundarySurface(BoundaryFace boundarySurface)
///	{// This function is probabbly not needed
///		// if the garbage collector takes care
///		// of loose pointers, but just in case ...
///		if(boundarySurface==null) return;
///		BoundaryFace	current=boundarySurface;
///		while (current.next!=null)
///		{	BoundaryFace	tmp=current;
///			current=current.next;
///			tmp.next=null;
///		}
///		boundarySurface=null;
///	}
///	public void setBoundarySurface
///	(	int[]	dim,
///		byte[]	grid,
///		double[]	X0,
///		double[]	DX
///	)
///	{	int 
///			n0=dim[0],
///			n1=dim[1],
///			n2=dim[2];
///		int[]	n=	new int[6];
///		double[]	x0 = new double[DIM];
///		double[] dx = new double[DIM];
///		for (int i=0;i<DIM;i++)
///		{	dx[i]=DX[i]/(double)dim[i];
///			x0[i]=X0[i]+0.5*dx[i];
///		}
///		m0=dim[0];m01=dim[0]*dim[1];
///		float[] v = new float[DIM];// face-nomral vector
///		// Delete the boundarySurface if not empty
///		if(boundarySurface!=null) 
///		{	deleteBoundarySurface(boundarySurface);
///			boundarySurface=null;
///		}
///		// Construct the boundarySurface
///		for (int i0=0;i0<n0;i0++)
///		for (int i1=0;i1<n1;i1++)
///		for (int i2=0;i2<n2;i2++)
///		{	int
///				nf, //number of boundarySurface faces of a cell
///				type=grid[I(i0,i1,i2)];
///			if(this.type!=type) continue;
///			if(i0>0   ) n[0]=I(i0-1,i1  ,i2  ); else n[0]=-1;
///			if(i0<n0-1) n[1]=I(i0+1,i1  ,i2  ); else n[1]=-1;
///			if(i1>0   ) n[2]=I(i0  ,i1-1,i2  ); else n[2]=-1;
///			if(i1<n1-1) n[3]=I(i0  ,i1+1,i2  ); else n[3]=-1;
///			if(i2>0   ) n[4]=I(i0  ,i1  ,i2-1); else n[4]=-1;
///			if(i2<n2-1) n[5]=I(i0  ,i1  ,i2+1); else n[5]=-1;
///			nf=0;
///			for (int i=0;i<DIM;i++) v[i]=0.0e0f;
///			for (int i=0;i<DIM;i++)
///			{	int	m=0;
///				for (int j=0;j<2;j++)
///				{	int	k=2*i+j;
///					if (n[k]<0||grid[n[k]]!=type)
///					{	// This is a boundarySurface cell:
///						// boundarySurface normal vector
///						// Add another boundarySurface element
///						// to the list
///						v[i]=-1.0f+2.0f*(float)j;
///						nf++;
///						m++;
///					}
///				}
///				if(m==2)
///				{	// If the boundarySurface was on the opposite 
///					// sides of the cell consider only one 
///					// in a positive direction.
///					v[i]=1.0f;
///				}
///			}
///			if (nf>0)
///			{	// Normalize to 1
///				double	a=0.0;
///				for (int i=0;i<DIM;i++)
///				{	double	r=v[i];
///					a+=r*r;
///				}
///				a=1./Math.sqrt(a);
///				for (int i=0;i<DIM;i++)
///					v[i]*=(float)a;
///				// Add the boundarySurface cell 
///				BoundaryFace newface = new BoundaryFace
///				(	(float)(x0[0]+dx[0]*i0),
///					(float)(x0[1]+dx[1]*i1),
///					(float)(x0[2]+dx[2]*i2),
///					v[0],v[1],v[2]
///				);
///				if (boundarySurface==null)
///					boundarySurface = newface;
///				else
///				{
///					newface.setNext(boundarySurface);
///					boundarySurface=newface;
///				}
///			}
///		}
///	}
///	public void	addBoundarySurfaceElement(BoundaryFace face)
///	{
///		if(boundarySurface==null)
///			boundarySurface = new BoundaryFace();
///	}
///	public void	printBoundarySurface()
///	{
///		for 
///		(	BoundaryFace element = boundarySurface;
///			element != null;
///			element=element.next
///		)
///		{	float[]
///				x=element.x,
///				v=element.v;
///			System.out.println
///			(	"X="+x[0]+","+x[1]+","+x[2]+"; V="+v[0]+","+v[1]+","+v[2]
///			);
///		}
///	}
///	public void saveMesh()
///	{
///	}
///}
