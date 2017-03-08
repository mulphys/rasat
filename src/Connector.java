
class Connector 
{	// Connects two components together
	int[][] ivars;
	Components[] comps;
	public Connector
	(
	 Components C1,
	 Components C2,
	 int[] ivar1,
	 int[] ivar2
	)
	{	comps = new Components[2];
		comps[0] = C1;
		comps[1] = C2;
		ivars = new int[2][];
		ivars[0] = new int[ivar1.length];
		ivars[1] = new int[ivar2.length];
		for(int i=0; i<ivar1.length; i++)
			ivars[0][i]=ivar1[i];
		for(int i=0; i<ivar2.length; i++)
			ivars[1][i]=ivar2[i];

	}
}
