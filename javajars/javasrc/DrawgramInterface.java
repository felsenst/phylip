package phylip;

import utilities.TestFileNames;
import drawgram.DrawgramUserInterface.DrawgramData;

import com.sun.jna.Library;
import com.sun.jna.Native;

public class DrawgramInterface {
    public interface Drawgram extends Library {
        public void  drawgram(
        		String  intree,
        		String  usefont,
        		String  plotfile,
        		String  plotfileopt,
        		String  treegrows,
        		String  treestyle,
        		boolean usebranchlengths,
        		double  labelangle,
        		boolean scalebranchlength,
        		double  branchlength,
        		double  breadthdepthratio,
        		double  stemltreedratio,
        		double  chhttipspratio,
        		double  xmarginratio,
        		double  ymarginratio,
        		String  ancnodes,
        		boolean doplot,
        		String  finalplotkind);

    }
	
	public boolean DrawgramRun(DrawgramData inVals){
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
			Drawgram Drawgram = (Drawgram) Native.loadLibrary("drawgram", Drawgram.class);
	        Drawgram.drawgram(
	        		inVals.intree,
	        		inVals.usefont,
	        		inVals.plotfile,
	        		inVals.plotfileopt,
	        		inVals.treegrows,
	        		inVals.treestyle,
	        		inVals.usebranchlengths,
	        		inVals.labelangle,
	        		inVals.scalebranchlength,
	        		inVals.branchlength,
	        		inVals.breadthdepthratio,
	        		inVals.stemltreedratio,
	        		inVals.chhttipspratio,
	        		inVals.xmarginratio,
	        		inVals.ymarginratio,
	        		inVals.ancnodes,
	        		inVals.doplot,
	        		inVals.finalplottype);
	        
			return true;
		}
    	catch(UnsatisfiedLinkError e)
    	{
     		new TestFileNames().LibraryMissing("DrawGram");
    		return false;
    	}
	}
}

	
