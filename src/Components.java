import java.awt.*;
import java.util.*;

public class Components
{
	static int	dim;
	public String name;
	public Color color;
	public BoundaryCell	boundaryCell=null; // Boundary vertex is a center
                                           // of the boundary cell
	public BoundaryContour[]	boundaryContours=null;

	public String toString()
	{
		return name;
	}
	// Handle boundary cells
	public void deleteBoundaryCells()
	{// This function is probabbly not needed
		// if the garbage collector takes care
		// of loose pointers, but just in case ...
		if(boundaryCell==null) return;
		BoundaryCell	current=boundaryCell;
		while (current.next!=null)
		{	BoundaryCell	tmp=current;
			current=current.next;
			tmp.next=null;
		}
		boundaryCell=null;
	}
	public void	addCell(BoundaryCell newCell)
	{
		if (boundaryCell==null)
			boundaryCell = newCell;
		else
		{	newCell.setNext(boundaryCell);
			boundaryCell=newCell;
		}
	}
	public void	printBoundaryCells()
	{
		BoundaryCell	curr=boundaryCell;
		for 
		(	curr=boundaryCell;
			curr!=null;
			curr=curr.next
		)
		{	int[] ind = curr.getInd();
			System.out.println("BoundaryCell: "+ind[0]+','+ind[1]+','+ind[2]);
		}
	}
	//Handle boundary contours
	public static void initBoundaryContourDim(int d)
	{	dim=d;
	}
	public void initBoundaryContours()
	{// This function is probabbly not needed
		// if the garbage collector takes care
		// of loose pointers, but just in case ...
		if(boundaryContours!=null)
		for (int i=0;i<dim;i++)
		{	if(boundaryContours[i]==null) continue;
			BoundaryContour	current=boundaryContours[i];
			while (current!=null)
			{	BoundaryContour	tmp=current;
				current=current.next;
				tmp.next=null;
			}
			boundaryContours[i]=null;
		}
		else
		{	boundaryContours = new BoundaryContour[dim];
			for (int i=0;i<dim;i++) boundaryContours[i]=null;
		}
	}
	public void initBoundaryContours(int iplane)
	{// This function is probabbly not needed
		// if the garbage collector takes care
		// of loose pointers, but just in case ...
		if(boundaryContours==null) 
		{	boundaryContours = new BoundaryContour[dim];
			for (int i=0;i<dim;i++)
				boundaryContours[i]=null;
		}
		if(boundaryContours[iplane]==null) return;
		BoundaryContour	current=boundaryContours[iplane];
		while (current.next!=null)
		{	BoundaryContour	tmp=current;
			current=current.next;
			tmp.next=null;
		}
		boundaryContours[iplane]=null;
	}
  	public void	addBoundaryContour
	(	int idir,// direction of the plane 
		int ipos,// position of the plane
		int nx,  // x-dimension of the plane
		int ny,  // y-dimension of the plane
		int ix0, //initial vertex x-position 
		int iy0, //initial vertex y-position
		byte[][][]	plane  
	)
	{	if (idir<0||idir>=dim)
		{	System.out.println("addContour: Error: plane orientation ="+idir);
			return;
		}
		int[]//neighbor cells indexes
			nbx= {-1,-1, 0, 1, 1, 1, 0,-1},
			nby= { 0, 1, 1, 1, 0,-1,-1,-1};
		int
			type=(int)plane[ix0][iy0][0],///PPP
			inb1=0,
			inb=inb1,
			nv=0, //number of vertexes in the contour
			ix1=ix0, iy1=iy0;
		LinkedList<Integer>	verts = new LinkedList<Integer>();
		((java.util.List<Integer>)verts).add(new Integer(((ix0-1)<<16)|(iy0-1)));nv++;///PPP
		plane[ix0][iy0][1]=1;//set status to first visited///PPP
		while(true)
		{// Construct the contour
			int ix=ix1, iy=iy1,jnb=0;
			//find the starting point
			for(;jnb<8;jnb++)
			{	inb=(inb1+jnb)%8;
				if(type!=(int)plane[ix1+nbx[inb]][iy1+nby[inb]][0])break;
			}
			if (jnb==8) return; // internal point
			for (jnb++;jnb<9;jnb++)
			{	inb=(inb1+jnb)%8;
				ix=ix1+nbx[inb];
				iy=iy1+nby[inb];
				if(type==(int)plane[ix][iy][0])break;
			}
			if(jnb==9) return; // only a single point
			if(ix==ix0&&iy==iy0)//visited: close the loop
			{	//Add vertex, and close the loop
				((java.util.List<Integer>)verts).add(new Integer(((ix-1)<<16)|(iy-1)));nv++;
				plane[ix][iy][1]++;///PPP
				if(plane[ix1][iy1][1]==1||plane[ix][iy][1]>1) break;
			}
			((java.util.List<Integer>)verts).add(new Integer(((ix-1)<<16)|(iy-1)));nv++;
			plane[ix][iy][1]++;//set status to visited///PPP
			ix1=ix; iy1=iy;
			inb1=(inb+4)%8;
		}
		// Add contour to the list
		if (boundaryContours[idir]==null)
			boundaryContours[idir] = new BoundaryContour(ipos,nv);
		else
		{	BoundaryContour	contour = new BoundaryContour(ipos,nv);
		///	contour.setNext(boundaryContours[idir]);
			contour.next=boundaryContours[idir];
			boundaryContours[idir]=contour;
		}
		ListIterator it = ((java.util.List)verts).listIterator();
		for (int i=0;i<nv;i++)
			boundaryContours[idir].putVertex(i,((Integer)it.next()).intValue());
	}
//-	public void	addBoundaryContour
//-	(	int idir,// direction of the plane 
//-		int ipos,// position of the plane
//-		int nx,  // x-dimension of the plane
//-		int ny,  // y-dimension of the plane
//-		int ix0, //initial vertex x-position 
//-		int iy0, //initial vertex y-position
//-		int type, // type of the component
//-		byte[][][]	plane  
//-	)
//-	{	if (idir<0||idir>=dim)
//-		{	System.out.println("addContour: Error: plane orientation ="+idir);
//-			return;
//-		}
//-		int[]//neighbor cells indexes
//-			nbx= {-1,-1, 0, 1, 1, 1, 0,-1},
//-			nby= { 0, 1, 1, 1, 0,-1,-1,-1};
//-		int
//-			inb0=0,
//-			inb=inb0,
//-			nv=0; //number of vertexes in the contour
//-		Object	verts = new LinkedList();
//-		while(true)
//-		{// Construct the contour
//-			int ix=ix0, iy=iy0,jnb=0;
//-			//find the starting point
//-			for(;jnb<8;jnb++)
//-			{	inb=(inb0+jnb)%8;
//-				if(type!=(int)plane[ix0+nbx[inb]][iy0+nby[inb]][0])break;
//-			}
//-			if (jnb==8) return; // internal point
//-			for (;jnb<8;jnb++)
//-			{	inb=(inb0+jnb)%8;
//-				ix=ix0+nbx[inb];
//-				iy=iy0+nby[inb];
//-				if(type==(int)plane[ix][iy][0])break;
//-			}
//-			if(jnb==8) return; // only a single point
//-			int status=(int)plane[ix][iy][1];
//-			if(status==1)//visited: close the loop
//-			{	//Add vertex
//-				((java.util.List)verts).add(new Integer(((ix-1)<<16)|(iy-1)));nv++;
//-				break;
//-			}
//-			if
//-			(	type!=(int)plane[ix+1][iy][0] || 
//-				type!=(int)plane[ix-1][iy][0] || 
//-				type!=(int)plane[ix][iy+1][0] || 
//-				type!=(int)plane[ix][iy-1][0]
//-			)
//-				((java.util.List)verts).add(new Integer(((ix-1)<<16)|(iy-1)));nv++;
//-			plane[ix][iy][1]=1;
//-			ix0=ix; iy0=iy;
//-			inb0=(inb+4)%8;
//-		}
//-		// Add contour to the list
//-		if (boundaryContours[idir]==null)
//-			boundaryContours[idir] = new BoundaryContour(ipos,nv);
//-		else
//-		{	BoundaryContour	contour = new BoundaryContour(ipos,nv);
//-		///	contour.setNext(boundaryContours[idir]);
//-			contour.next=boundaryContours[idir];
//-			boundaryContours[idir]=contour;
//-		}
//-		ListIterator it = ((java.util.List)verts).listIterator();
//-		for (int i=0;i<nv;i++)
//-			boundaryContours[idir].putVertex(i,((Integer)it.next()).intValue());
//-	}
	public void	printBoundaryContours()
	{	if (boundaryContours==null) return;
		for (int i=0;i<dim;i++)
		{	BoundaryContour	curr;
			System.out.println("BoundaryContours in plane "+i);
			for 
			(	curr=boundaryContours[i];
				curr!=null;
				curr=curr.next
			)
			{	int n = curr.Length();
				System.out.println("\t Contour length: "+n);
				for (int j=0;j<n;j++)
					System.out.println("\t\t"+curr.ind[j]);
			}
		}
	}
}
