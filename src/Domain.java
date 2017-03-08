import java.io.*;
import java.util.Vector;
import java.net.URL;

public class	Domain 
{
	boolean	modified=false;
	static int
		DIM,
		m0,m01; // used in I(i,j,k)
	int[]
		Dim, // dimensions of the grid
		Iorigin; // Offset from the space origin
	double[]	X0,DX;//coordinates origin,
	               // and domain dimensions
	public String	name;
	String path=".";
	public byte[]	grid=null;
	DomainList	neibs; // domain neighbors
	public Domain	next,prev;

	static public int	I
	(	int i0, int i1, int i2
	)
	{	return	i0+i1*m0+i2*m01;
	}
	public Domain
	(	int	i0,
		int	i1,
		int	i2,
		Space space
	)
	{	DIM = space.getDim();
		Dim = new int[DIM];
		Iorigin = new int[DIM];
		Iorigin[0]=i0;
		Iorigin[1]=i1;
		Iorigin[2]=i2;
		X0 = new double[DIM];
		DX = new double[DIM];
		for(int i=0;i<DIM;i++)
		{	X0[i]=space.getX0(i)+(double)Iorigin[i]*space.getdX(i);
			DX[i]=0.0;
		}
		name=null;
		grid=null;
		neibs=null;
		next=prev=this;
		modified=false;
		path=space.getPath();
	}
	public	Domain 
	(
		String name,
		int i0, int i1, int i2,
		Space space
	)
	{
		this(i0, i1, i2, space);
		this.name=name;
	}
	public	Domain 
	(
		String name,
		int i0, int i1, int i2,
		int n0, int n1, int n2,
		Space space
	)
	{	this(name,i0,i1,i2,space);
		Dim[0]=n0;
		Dim[1]=n1;
		Dim[2]=n2;
		m0=Dim[0];m01=Dim[0]*Dim[1];
		for (int i=0;i<DIM;i++)
			DX[i]=space.getdX(i)*Dim[i];
	}
	public void	LoadGhost (String url_str)
	{	if(name==null)
		{	Warning("Load: Domain name unknown");
			return;
		}
		if(grid!=null)
		{	Warning("Grid's not empty");
			return;
		}
		if(Dim==null)
		{	Dim=new int[DIM];
		}
		try
		{	//Open domain file and read 
			URL file_url = new URL(url_str+name+".dom");
			InputStream is = file_url.openStream();
			DataInputStream	inp = new DataInputStream
			(	new BufferedInputStream(is)
			);
			for (int i=0;i<DIM;i++) Dim[i]=inp.readInt();
			for (int i=0;i<DIM;i++) X0[i]=inp.readDouble();
			for (int i=0;i<DIM;i++)	DX[i]=inp.readDouble();
			inp.close();
		}	catch(Exception ex)
		{	Warning(ex.toString());
		}
		modified=false;
	}
	public void	Load (String url_str)
	{	if(name==null)
		{	Warning("Load: Domain name unknown");
			return;
		}
		if(grid!=null)
		{	Warning("Grid's not empty");
			return;
		}
		if(Dim==null) Dim=new int[DIM];
		if(X0==null) X0=new double[DIM];
		if(DX==null) DX=new double[DIM];
		if(modified==true)
		{	Warning("Domain changed, save before loading");
			return;
		}
		try
		{
			URL file_url = new URL(url_str+name+".dom");
			InputStream is = file_url.openStream();
			DataInputStream	inp = new DataInputStream
			(	new BufferedInputStream(is)
			);
			for (int i=0;i<DIM;i++) Dim[i]=inp.readInt();
			Activate();
			for (int i=0;i<DIM;i++) X0[i]=inp.readDouble();
			for (int i=0;i<DIM;i++)	DX[i]=inp.readDouble();
			for (int i0=0;i0<Dim[0];i0++)
			for (int i1=0;i1<Dim[1];i1++)
			for (int i2=0;i2<Dim[2];i2++)
				grid[I(i0,i1,i2)]=inp.readByte();
			inp.close();
		}	catch(Exception ex)
		{	Warning(ex.toString());
		}
		modified=false;
	}
	public void error(String message)
	{
		System.out.println(message);
	}
	public void Exit ()
	{
///		Save();
		Deactivate();
	}
	public void	setName
	(	String name
	)
	{	this.name=name;
	}
	public void	setModified ()
	{	modified=true;
	}
	public void	setType
	(
		int i, int j, int k,
	        byte type
	)
	{	if (grid==null) 
		{	Warning("Can't set type to an empty grid");
			return;
		}
		try
		{	///grid[i][j][k].setType(type);
			grid[I(i,j,k)]=type;
		}	catch (Exception e)
		{
			Warning("Can't set type at cell "+i+","+j+","+k);
			return;
		}
	}
	public void setBallType
	(
		int i, int j, int k, // center coordinats
		int	r,	  //radius
		byte	type
	)
	{	if (grid==null) 
		{	Warning("Can't set type to an empty grid");
			return;
		}
		try
		{	int	ni=Dim[0],nj=Dim[1],nk=Dim[2];
			for (int ii=i-r;ii<i+r;ii++)
			{	if(ii<0||ii>ni) continue;
			for (int jj=j-r;jj<j+r;jj++)
			{	if(jj<0||jj>nj) continue;
			for (int kk=k-r;kk<k+r;kk++)
			{	if(kk<0||kk>nk) continue;
				float	
					x=(float)(ii-i),
					y=(float)(jj-j),
					z=(float)(kk-k),
					d=(float)Math.sqrt(x*x+y*y+z*z);
				if(d>(float)r) continue;
				grid[I(ii,jj,kk)]=(byte)type;
			}
			}
			}
		}	catch (Exception e)
		{
			Warning("Can't set ball type at cell "+i+","+j+","+k);
			return;
		}
	}
	public byte	getType
	(
		int i, int j, int k
	)
	{	byte	type=0;
		if (grid==null) 
		{	Warning("Can't set type to an empty grid");
			return 0;
		}
		try
		{	type=grid[I(i,j,k)];
		}	catch (Exception e)
		{	Warning("Can't set type at cell "+i+","+j+","+k);
			return 0;
		}
		return type;
	}
	public void	Activate ()
	{
		if(Dim==null)return;
		if(grid!=null)return;
		m0=Dim[0];m01=Dim[0]*Dim[1];
		try
		{	grid = new byte[Dim[0]*Dim[1]*Dim[2]];
			for (int i=0;i<Dim[0];i++)
			for (int j=0;j<Dim[1];j++)
			for (int k=0;k<Dim[2];k++)
				grid[I(i,j,k)]=0;
		}
		catch (Exception e)
		{	grid=null;
			Warning("Can't activate domain"+name);
		}
	}
	public void	Deactivate()
	{	if (modified==true) 
		{	Warning("Domain is modified: use Save or Delete first");
			return;
		}
		grid=null;
	}
	public void	deleteGrid()
	{
		grid=null;
		Dim[0]=Dim[1]=Dim[2]=0;
		m0=Dim[0];m01=Dim[0]*Dim[1];
	}
	public void Print()
	{
		System.out.println("Domain name = " + name);
		if(Dim!=null)
			System.out.println("\tDimensions: "+Dim[0]+","+Dim[1]+","+Dim[2]);
		if(grid!=null)
			for(int i=0;i<getDim(0);i++)
			for(int j=0;j<getDim(1);j++)
			for(int k=0;k<getDim(2);k++)
				System.out.println("\t"+i+":"+j+":"+k+": "+grid[I(i,j,k)]);
		else
			System.out.println("\tGhost domain");
	}
	public String getName()
	{
		return name;
	}
	public boolean isoverlap(Domain a)
	{
		int[][] ia = new int[2][DIM];
		int[][] ib = new int[2][DIM];
		int[][] iv = new int[2*2*2][DIM];
		Domain b=this;
		do
		{	for (int i=0;i<DIM;i++)
			{
				ia[0][i]=a.getIorigin(i);
				ib[0][i]=b.getIorigin(i);
				ia[1][i]=ia[0][i]+a.getDim(i);
				ib[1][i]=ib[0][i]+b.getDim(i);
			}
			for (int n=0;n<DIM;n++)	
			for (int i=0;i<2;i++)
			for (int j=0;j<2;j++)
			for (int k=0;k<2;k++)
			{	int	m=0;
				if(n==0&&i==1||n==1&&j==1||n==2&&k==1)m=1;
				iv[i+2*j+4*k][n]=ia[m][n];
			}
			for (int i=0;i<8;i++)
			{	int inside=0;
				for (int j=0;j<DIM;j++)
				if(iv[i][j]>=ib[0][j]&&iv[i][j]<=ib[1][j])
					inside++;
				if(inside==DIM)
					return true;
			}
			b=b.next;
		}	while(b!=this);
		return false;
	}
	public void Append(Domain newmember)
	{
		if(this.isoverlap(newmember)==true)
		{	Warning("Domain Append Failed: OVERLAP DETECTED");
			return;
		}
		Domain	tmp=prev;
		prev.next=newmember;
		prev=newmember;
		newmember.prev=tmp;
		newmember.next=this;
	}
	int[]	getDim()
	{
		return Dim;
	}
	int	getDim(int i)
	{
		return Dim[i];
	}
	int	getIorigin(int i)
	{
		return Iorigin[i];
	}
	byte[]	getGrid()
	{	return	grid;
	}
///	void	setPlane(int iplane)
///	{	this.iplane=iplane;
///	}
///	int	getPlane()
///	{	return iplane;
///	}
///	void	setPlanePos(int ipos)
///	{
///		if(ipos >= Dim[iplane]) iplanepos=Dim[iplane]-1;
///		else iplanepos=ipos;
///	}
///	int	getPlanePos()
///	{	return	iplanepos;
///	}
///	int	getPlaneNx()
///	{	int	n=0;
///		switch(iplane)
///		{	case 0: n=Dim[1]; break;
///			case 1: n=Dim[0]; break;
///			case 2: n=Dim[0]; break;
///		}
///		return n;
///	}
///	int	getPlaneNy()
///	{	int	n=0;
///		switch(iplane)
///		{	case 0: n=Dim[2]; break;
///			case 1: n=Dim[2]; break;
///			case 2: n=Dim[1]; break;
///		}
///		return n;
///	}
	void	Warning(String message)
	{	System.out.println(message);
	}
	void	Connect(Domain neibdom)
	{
		if (neibdom==this)return;
		// Check if already connected
		if(neibs!=null)
		{	DomainList	neib=neibs;
			do
			{	if(neib.dom==neibdom)
				{	Warning("Domains already connected");
					return;
				}
				neib=neib.next;
			}	while (neib!=neibs);
			neibs.Append(neibdom);
		}
		else
			neibs = new DomainList(neibdom);
	}
	void	Disconnect(Domain neibdom)
	{	// Remove neibdom from the neighbor list
		if (neibdom==this)return;
		// Check if already connected
		if(neibs==null)return;
		DomainList	neib=neibs;
		do
		{	if(neib.dom==neibdom)
			{	Warning("Disconnecting domain "+neib.dom.getName());
				if(neibs.next==neibs)
				{	neibs.next=neibs.prev=null;
					neibs.dom=null;
					neibs=null;
					return;
				}
				if(neib==neibs)neibs=neibs.next;
				neib.prev.next=neib.next;
				neib.next.prev=neib.prev;
				neib.dom=null;
				neib.next=neib.prev=null;
				neib=null;
				return;
			}
			neib=neib.next;
		}	while (neib!=neibs);
	}
	void	Disconnect()
	{	// Remove all neighbors
		if(neibs==null)return;
		if(neibs.next!=neibs)
		{	DomainList	neib=neibs.next;
			do
			{	neib.dom.Disconnect(this);
				neib.dom=null;
				DomainList old=neib;
				neib=neib.next;
				old.prev=old.next=null;
				old=null;
			}	while(neib!=neibs);
		}
		neibs.dom.Disconnect(this);
		neibs.dom=null;
		neibs.next=neibs.prev=null;
		neibs=null;
	}
	public byte[]	getgrid()
	{	return grid;
	}
	public double[]	getX0()
	{	return X0;
	}
	public double	getX0(int i)
	{	return X0[i];
	}
	public double[]	getDX()
	{	return DX;
	}
	public void setBoundaryCells
	(
		Vector	vComponents
	)
	{	int 
			n0=Dim[0],
			n1=Dim[1],
			n2=Dim[2],
			maxtypes=vComponents.size();
		int[]	n=	new int[6];
		m0=n0;m01=n0*n1;
		for(int i=0;i<vComponents.size();i++)
		{	// Delete boundaryCells if not empty
			((Components)vComponents.elementAt(i)).deleteBoundaryCells();
		}
		BoundaryCell.init(DIM,Dim,X0,DX);
		for (int i0=0;i0<n0;i0++)
		for (int i1=0;i1<n1;i1++)
		for (int i2=0;i2<n2;i2++)
		{	int	type=(int)grid[I(i0,i1,i2)];
			if(type==0)continue;
			if(i0>0   ) n[0]=I(i0-1,i1  ,i2  ); else n[0]=-1;
			if(i0<n0-1) n[1]=I(i0+1,i1  ,i2  ); else n[1]=-1;
			if(i1>0   ) n[2]=I(i0  ,i1-1,i2  ); else n[2]=-1;
			if(i1<n1-1) n[3]=I(i0  ,i1+1,i2  ); else n[3]=-1;
			if(i2>0   ) n[4]=I(i0  ,i1  ,i2-1); else n[4]=-1;
			if(i2<n2-1) n[5]=I(i0  ,i1  ,i2+1); else n[5]=-1;
			for (int i=0;i<DIM;i++)
			{	int	j=0;
				for (;j<2;j++)
				{	int	k=2*i+j;
					if (n[k]<0||grid[n[k]]!=type)
					{	// This is a boundary cell:
						// Add it to the list of boundary vertexes
						if((int)type>maxtypes)
						{	error("Unregistered component type="+type+" in cell "+i0+','+i1+','+i2);
							continue;
						}
						((Components)vComponents.elementAt(type)).addCell
						(	new BoundaryCell(i0,i1,i2)
						);
						break;
					}
				}
				if(j<2) break;
			}
		}
	}
	public void setBoundaryContours
	(
		int iplane, // cross-sectional plane orientation
		int dz, // plane position increment
		Vector	vComponents
	)
	{	int 
			nx,ny,nz,
			n0=Dim[0],
			n1=Dim[1],
			n2=Dim[2],
			maxtypes=vComponents.size();
		m0=n0;m01=n0*n1;
		Components.initBoundaryContourDim(DIM);
		// Initialize class:
///		BoundaryContour.init(DIM,Dim,X0,DX);
		for(int i=0;i<vComponents.size();i++)
		{	// Delete boundaryContours if not empty
			((Components)vComponents.elementAt(i)).initBoundaryContours(iplane);
		}
		// Select a plane
		switch(iplane)
		{	case 0:
			nx=n1; ny=n2; nz=n0;
			break;
			case 1:
			nx=n2; ny=n0; nz=n1;
			break;
			case 2:
			nx=n0; ny=n1; nz=n2;
			break;
			default:
			    //SystemLog.displayMessage("setBoundaryContours: Error: iplane="+iplane);
			return;
		}
		byte[][][]	plane = new byte[nx+2][ny+2][2];
		for (int i=0;i<nx+2;i++)
		{	plane[i][0][0]=0;
			plane[i][ny+1][0]=0;
		}
		for (int i=0;i<ny+2;i++)
		{	plane[0][i][0]=0;
			plane[nx+1][i][0]=0;
		}
		for (int iz=0;iz<nz;iz+=dz)
		{	int	i0,i1,i2;
			switch(iplane)
			{	case 0:
					i0=iz;
					for (i1=0;i1<n1;i1++)
					for (i2=0;i2<n2;i2++)
					{	plane[i1+1][i2+1][0]=grid[I(i0,i1,i2)];
						plane[i1+1][i2+1][1]=0;
					}
				break;
				case 1:
					i1=iz;
					for (i0=0;i0<n0;i0++)
					for (i2=0;i2<n2;i2++)
					{	plane[i2+1][i0+1][0]=grid[I(i0,i1,i2)];
						plane[i2+1][i0+1][1]=0;
					}
				break;
				case 2:
					i2=iz;
					for (i0=0;i0<n0;i0++)
					for (i1=0;i1<n1;i1++)
					{	plane[i0+1][i1+1][0]=grid[I(i0,i1,i2)];
						plane[i0+1][i1+1][1]=0;
					}
				break;
				default:
				    //SystemLog.displayMessage("setBoundaryContours: Error: iplane="+iplane);
				return;
			}
			setBoundaryContours
			(	iplane, iz, nx, ny, plane,
				vComponents
			);
		}
	}
	public void setBoundaryContours
	(	int iplane, int iz, // plane orientation and position
		int nx, int ny, // plane dimenstions
		byte[][][]	plane,
		Vector	vComponents
	)
	{	int	max_types=vComponents.size();
		for(int iy=1;iy<=ny;iy++)
		for(int ix=1;ix<=nx;ix++)
		{	int
				type=(int)plane[ix][iy][0],
				status=(int)plane[ix][iy][1];
			if(type<=0)continue;
			if(type>=max_types)
			{	
				Warning("setContours: Wrong type detected " + type + " >= " + max_types);
				continue;
			}
			///if(plane[ix-1][iy][0]==type||status==1) continue;
			if(plane[ix-1][iy][0]==type||status!=0) continue;
			((Components)vComponents.elementAt(type)).addBoundaryContour
			(iplane, iz, nx, ny, ix, iy, plane);
		}
	}

}
class	DomainList
{
	Domain	dom;
	DomainList	next,prev;
	public DomainList(Domain dom)
	{
		this.dom=dom;
		next=prev=this;
	}
	public void Append(Domain newdom)
	{
		DomainList
			newmember = new DomainList(newdom),
			tmp=prev;
		prev.next=newmember;
		prev=newmember;
		newmember.prev=tmp;
		newmember.next=this;
	}
}

