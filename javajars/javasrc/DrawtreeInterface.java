package phylip;

import utilities.TestFileNames;
import drawtree.DrawtreeUserInterface.DrawtreeData;

import com.sun.jna.Library;
import com.sun.jna.Native;

public class DrawtreeInterface {
    public interface Drawtree extends Library {
        public void  drawtree(
        	 	String  intree,
        		String  plotfile,
          		String  plotfileopt,
          	  	String  usefont,
        		String  treegrows,
        		boolean usebranchlengths,
        		String  labeldirec,
        		Double  labelangle,
        		Double  treerotation,
        		Double  treearc,
        		String  iterationkind,
        		int     iterationcount,
        		boolean regularizeangles,
        		boolean avoidlabeloverlap,
        		boolean branchrescale,
        		Double  branchscaler,
        		Double  relcharhgt,
        		double  xmarginratio,
        		double  ymarginratio,
        		boolean doplot,
        		String  finalplotkind);
    }
	
	public boolean DrawtreeRun(DrawtreeData inVals){
		TestFileNames test = new TestFileNames();
		
		if (!test.FileAvailable(inVals.intree, "Intree"))
		{
			return false;
		}
		
		if (inVals.doplot) // only check if final plot
		{ 
			String opt = test.FileAlreadyExists(inVals.plotfile, "Plotfile");
			if (opt == "q")
			{
				return false;
			}
			else
			{
				inVals.plotfileopt = opt;
			}
		}
		
		// at this point we hook into the C code
		try
		{		
			Drawtree Drawtree = (Drawtree) Native.loadLibrary("drawtree", Drawtree.class);
	        Drawtree.drawtree(
	        		inVals.intree,
	        		inVals.plotfile,
	        		inVals.plotfileopt,
	        		inVals.usefont,
	        		inVals.treegrows,
	        		inVals.usebranchlengths,
	        		inVals.labeldirec,
	        		inVals.labelangle,
	        		inVals.treerotation,
	        		inVals.treearc,
	        		inVals.iterationkind,
	        		inVals.iterationcount,
	        		inVals.regularizeangles,
	        		inVals.avoidlabeloverlap,
	        		inVals.branchrescale,
	        		inVals.branchscaler,
	        		inVals.relcharhgt,
	        		inVals.xmarginratio,
	        		inVals.ymarginratio,
	        		inVals.doplot,
	        		inVals.finalplottype);
	        return true;
		}
		catch(UnsatisfiedLinkError e)
		{
	 		new TestFileNames().LibraryMissing("DrawTree");
			return false;
		}
	}
}
