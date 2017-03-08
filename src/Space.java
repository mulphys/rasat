import java.io.*;
import java.util.Vector;

public class Space
{
	final static int DIM = 3;
	double[]	X0,dX;//coordinates origin,
	               // and cell-size dimensions
	public Vector<Components> vComponents; //Store the component IDs

	static String	Name,Path;
	static Connector[]	connectors;
///	DomainComponent	comps; // collection of components
	public static Domain	first_domain;

	public	Space
	(	String name,
		double	d0,
		double	d1,
		double	d2
	)
	{
		Name=name;
		X0 = new double[DIM];
		dX = new double[DIM];
		for (int i=0;i<DIM;i++)
		{
			X0[i]=0.0; // default origin
		}
		dX[0]=d0;
		dX[1]=d1;
		dX[2]=d2;
		Path=".";
		vComponents = new Vector<Components>(); //Store the component IDs
	}
	public String getName()
	{
		return Name;
	}
	public String getPath()
	{
		return Path;
	}
	public void setPath(String path)
	{
		Path=path;
	}
	public int getDim()
	{
		return DIM;
	}
	public double getX0(int i)
	{
		return X0[i];
	}
	public double getdX(int i)
	{
		return dX[i];
	}
	public double[] getdX()
	{
		return dX;
	}
	public Domain getHead()
	{
		return first_domain;
	}
	public Vector getComponents()
	{
		return vComponents;
	}
	public static void main (String[] argv)
	{	int	ndom=3,
			nx=2,ny=3,nz=4,
			ix=0,iy=0,iz=0;
		Space	space = new Space("Test",1.0,1.0,1.0);
		System.out.println("Creating domains:\n");
		// Initializing first domain
		first_domain = new Domain("Domain0",ix,iy,iz,nx,ny,nz,space);
		first_domain.Print();
		// Initializing more domains
		for (int i=1;i<ndom;i++)
		{	nx+=i;
			ny+=i;
			nz+=i;
			ix+=nx;
			System.out.println("Creating domain:\n" + i);
			Domain new_domain = new Domain("Domain"+i,ix,iy,iz,nx,ny,nz,space);
			first_domain.Append(new_domain);
		}
		System.out.println("Output domains:\n");
		Domain	
			head=space.getHead(),
			dom=head;
		do
		{	dom.Print();
			dom=dom.next;
		}	while (dom!=head);
		// Create grids on each domain
		dom=first_domain;
		do
		{	int
				ni=dom.getDim(0),
				nj=dom.getDim(1),
				nk=dom.getDim(2);
			dom.Activate();
			// Fill in cubicles:
			for (int i=0;i<ni;i++)
			for (int j=0;j<nj;j++)
			for (int k=0;k<nk;k++)
			{	byte type=1; // assign types (colors) to the grid
				if(k<nk/2)type=2;
				dom.setType(i,j,k,type);
			}
			dom.setModified();
			System.out.println("Printing domain "+dom.name);
			dom.Print();
			dom.Exit();// Domain remains in the list but
			           // grid is dumped into the file
			           // and memory released
			dom=dom.next;
		}	while (dom!=first_domain);
		// Load one domain at a time and display
		System.out.println("Loading and printing domains");
		dom=first_domain;
		do
		{	dom.Load("");
			dom.Print();
			dom=dom.next;
		}	while(dom!=first_domain);
		// Connecting all domains
		dom=first_domain;
		do
		{	Domain nextdom=dom;
			do
			{	Connect(dom,nextdom);
				nextdom=nextdom.next;
			}	while(nextdom!=dom);
			dom=dom.next;
		}	while(dom!=first_domain);
		// Print connectivity
		System.out.println("Domain connectivity");
		dom=first_domain;
		do
		{	System.out.println(dom.getName());
			DomainList neibs=dom.neibs,neib=neibs;
			if(neibs!=null)
			do
			{	System.out.println("\t"+neib.dom.getName());
				neib=neib.next;
			}	while(neib!=neibs);
			dom=dom.next;
		}	while(dom!=first_domain);
		// Save domain connectivity
		saveCon();
		// Create components
///		for(int i=0;i<3;i++)
///		{	comp = new DomainComponentID();
///			comp.name="Comp"+i;
///			comp.color=new Color(10,100,200);
///			// Create component boundary
///			comp.setBoundarySurface
///			(
///				first_domain.getDim(),
///				first_domain.getgrid(),
///				first_domain.getX0(),
///				first_domain.getDX()
///			);
///			comp.printBoundarySurface();
///		}
		// Deleting domains
		System.out.println("Deleting domains");
		for (int i=0;i<ndom;i++)
		{	Remove(first_domain);
			// Print domains
			System.out.println("Remaining domains:");
			if(first_domain!=null)
			{	dom=first_domain;
				do
				{	System.out.println(dom.getName());
					DomainList neibs=dom.neibs,neib=neibs;
					if(neibs!=null)
					do
					{	System.out.println("\t"+neib.dom.getName());
						neib=neib.next;
					}	while(neib!=neibs);
					dom=dom.next;
				}	while(dom!=first_domain);
			}
		}
	}

	public static void	RemoveDomains()
	{
		while(first_domain!=null) Remove(first_domain);
	}
	static void	Remove(Domain dom)
	{
		if (first_domain==null) return;
		if (first_domain.next==first_domain)
		{
			first_domain.neibs=null;
			first_domain.next=first_domain.prev=null;
			first_domain=null;
			return;
		}
		if (first_domain==dom) first_domain=first_domain.prev;
		dom.Disconnect();
		dom.prev.next=dom.next;
		dom.next.prev=dom.prev;
		dom.next=dom.prev=null;
		dom.X0=dom.DX=null;
		dom.grid=null;
		dom.Dim=null;
		dom=null;
	}
	static void	Connect(Domain D1, Domain D2)
	{
		D1.Connect(D2);
		D2.Connect(D1);
	}
	static void	Disconnect(Domain D1, Domain D2)
	{
		D1.Disconnect(D2);
		D2.Disconnect(D1);
	}
	static void	saveCon()
	{	try
		{	File	outFile = new File(Name+".con");
			FileOutputStream	outStream = new FileOutputStream(outFile);
			PrintWriter	out = new PrintWriter(outStream);
			Domain	dom=first_domain;
			do
			{	out.println(dom.getName());
				DomainList neibs=dom.neibs,neib=neibs;
				if(neibs!=null)
				do
				{	out.println("\t"+neib.dom.getName());
					neib=neib.next;
				}	while(neib!=neibs);
				dom=dom.next;
			}	while(dom!=first_domain);
			out.close();
		}	catch(Exception ex)
		{	System.out.println(ex.toString());
		}
	}
	static void	readCon() 
	{	try
		{	File	inpFile = new File(Name+".con");
			FileInputStream	inpStream = new FileInputStream(inpFile);
			BufferedInputStream inp = new BufferedInputStream(inpStream);
			if (first_domain==null)return;
			Domain	dom=first_domain;
			do
			{
///DDD: Could you finish this for me? Something is wrong in
/// 	the next line. If you can fix it the rest should be fine.
		//	inp.readLine(dom.getName());
		//		DomainList neibs=dom.neibs,neib=neibs;
		//		if(neibs!=null)
		//		do
		//		{	out.println("\t"+neib.dom.getName());
		//			neib=neib.next;
		//		}	while(neib!=neibs);
				dom=dom.next;
			}	while(dom!=first_domain);
			inp.close();
		}	catch(Exception ex)
		{	System.out.println(ex.toString());
		}
	}
}
