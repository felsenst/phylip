package phylip;
//import javax.swing.SwingUtilities;

/*
import java.awt.Color;
import java.awt.Dimension;
import java.awt.EventQueue;
import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;

import javax.swing.JButton;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JScrollPane;
import javax.swing.JTextField;
import javax.swing.SwingConstants;
*/
import org.eclipse.swt.*;
import org.eclipse.swt.widgets.*;

//import CliqueOldUserInterface.Clique;
//import utilities.DisplayProgress;
import utilities.TestFileNames;

import org.eclipse.swt.layout.*;
import org.eclipse.swt.events.*;
//import org.eclipse.swt.graphics.*;
//import org.eclipse.swt.internal.cocoa.NSFont;

import com.sun.jna.Library;
import com.sun.jna.Native;

//import CliqueOldUserInterface.Clique;

public class CliquePrototypeInterface {
	  public interface Clique extends Library {
	        public void  clique(
	        		String infile,
	        		String outfile,
	        		String outfileopt,
	        		String outtree,
	        		String outtreeopt,
	        		String factfile,
	        		String ancfile,
	        		String weightfile,
	         	boolean useAncStates,
	        		boolean useFactors,
	        		boolean useWeights,
	        		int numWeights,
	        		boolean noCliqueSize,
	        		int cliqueSize,
	        		boolean multipleDataSets,
	        		boolean multipleWeightSets,
	        		boolean useOutGroupRoot,
	        		int outRoot,
	        		boolean printRawData,
	        		boolean progressBar,
	        		boolean printMatrix,
	        		boolean printTree,
	        		boolean outputTree
	    	        );
	    	    }

	//public class cliqueData{
		String infile;
		String outfile;
		String outfileopt;
		String outtree;
		String outtreeopt;
		String factfile;
		String ancfile;
		String weightfile;
		String compMatrix;
		boolean useAncStates;
		boolean useFactors;
		boolean useWeights;
		int numWeights;
		boolean noCliqueSize;
		int cliqueSize;
		boolean multipleDataSets;
		boolean multipleWeightSets;
		boolean useOutGroupRoot;
		int outRoot;
		boolean printRawData;
		boolean progressBar;
		boolean printMatrix;
		boolean printTree;
		boolean outputTree;	
	//}
	
	//private static cliqueData clique;

	private static Display display;
	private static Shell shShell;
	private static Button btnInputFile;
	private static Text txtInputFile;
	
	private static Button btnOutputFile;
	private static Text txtOutputFile;
	
	private static Button rdbtnOtree;
	private static Text txtOuttree;
	private static Button btnOutTree;
	
	private static Button rdbtnFac;
	private static Text txtFactfile;
	
	private static Button rdbtnAnc;
	private static Text txtAncfile;
	
	private static Button rdbtnWgt;
	private static boolean weightstate;
	private static Text txtWgtsfile;
	
	private static Button rdbtnOgrp;
	private static Text txtOutRoot;
	private static Button rdbtnOutGroup;
	private static boolean outspecies;
	
	private static Button rdbtnCliq;
	private static Label lblSize; 
	private static Text txtCliqueSize;
	
	private static Button rdbtnBrun;
	private static Button rdbtnCmat;
	private static Button rdbtnPtree;
	
	private static Label lblNumSets;
	private static Text txtSets;
	
	private static Label lblMultipleDatasets;
	private static Button rdbtnMdsY;
	private static Button rdbtnMdsN;
	private static Button rdbtnMwsY;

	private static Button rbtnExecute;
	private static Button rbtnQuit;
	
	private static Button rdbtnPrintInd;
	
	// event handlers
	protected static void AncToggle() {	
		if (txtAncfile.getEnabled())
		{
			txtAncfile.setEnabled(false);
			txtAncfile.setText("ancfile");
		}
		else
		{
			txtAncfile.setEnabled(true);			
			FileDialog fd = new FileDialog(shShell, SWT.OPEN);
			String selected = fd.open();
			if (selected != null)
			{
				txtAncfile.setText(selected);
			}
		}
	}
	
	protected static void FacToggle() {		
		if (txtFactfile.getEnabled())
		{
			txtFactfile.setEnabled(false);
			txtFactfile.setText("factfile");
		}
		else
		{
			txtFactfile.setEnabled(true);			
			FileDialog fd = new FileDialog(shShell, SWT.OPEN);
			String selected = fd.open();
			if (selected != null)
			{
				txtFactfile.setText(selected);
			}
		}
	}

	protected static void WgtToggle() {		
		if (txtWgtsfile.getEnabled())
		{
			txtWgtsfile.setEnabled(false);
			txtWgtsfile.setText("weightfile");
			weightstate = false;
		}
		else
		{/*
			if (rdbtnMdsY.getEnabled()){
				rdbtnMwsY.setEnabled(false);
				MessageBox messageBox = new MessageBox(shShell, SWT.ICON_WARNING | SWT.OK );
		        
		        messageBox.setText("Warning");
				String msg = "Multiple data sets take precidence over multiple weights. \"Yes\" cannot be selected.";
		        messageBox.setMessage(msg);
		        //int buttonID = messageBox.open();
		        return;
	        }
	        */
			txtWgtsfile.setEnabled(true);			
			FileDialog fd = new FileDialog(shShell, SWT.OPEN);
			String selected = fd.open();
			if (selected != null)
			{
				txtWgtsfile.setText(selected);
			}
			weightstate = true;
		}
	}
	
	protected static void CliqSizeToggle() {		
		if (txtCliqueSize.getEnabled())
		{
			txtCliqueSize.setEnabled(false);
			txtCliqueSize.setText("1");
		}
		else
		{
			txtCliqueSize.setEnabled(true);	
		}
	}

	
	protected static void OutGroupToggle() {		
		if (rdbtnOutGroup.getEnabled()){
			rdbtnOutGroup.setEnabled(false);
			txtOutRoot.setEnabled(false);
			txtOutRoot.setText("1");
		}	
		else{
			rdbtnOutGroup.setEnabled(true);
			txtOutRoot.setEnabled(true);
		}		
	}
	
	protected static void SpeciesToggle() {	
		if (outspecies){
			rdbtnOutGroup.setText("Outgroup species:");
		}	
		else{
			rdbtnOutGroup.setText("At species number:");
		}		
		outspecies = !outspecies;
	}
	
	protected static void MultiSetToggle(int kind) {	
		
		if (kind == 1){
			// Wgttype.DATA
			rdbtnMdsY.setSelection(true);
			rdbtnMdsN.setSelection(false);
			rdbtnMwsY.setSelection(false);
			lblNumSets.setEnabled(true);
			txtSets.setEnabled(true);
			if (weightstate == true){
				rdbtnWgt.setSelection(true);
				txtWgtsfile.setEnabled(true);
			}			
			else{
				rdbtnWgt.setSelection(false);
				txtWgtsfile.setEnabled(false);
			}	
		}	
		else{
			if (kind == 0){
				// Wgttype.NONE
				rdbtnMdsY.setSelection(false);
				rdbtnMdsN.setSelection(true);
				rdbtnMwsY.setSelection(false);	
				lblNumSets.setEnabled(false);
				txtSets.setEnabled(false);
				
				if (weightstate == true){
					rdbtnWgt.setSelection(true);
					txtWgtsfile.setEnabled(true);
				}
				else{
					rdbtnWgt.setSelection(false);
					txtWgtsfile.setEnabled(false);					
				}
			}
			else {
				// (kind == Wgttype.WEIGHT == 2)
				rdbtnMdsY.setSelection(false);
				rdbtnMdsN.setSelection(false);
				rdbtnMwsY.setSelection(true);
				lblNumSets.setEnabled(true);
				txtSets.setEnabled(true);
				
				/*
				if (rdbtnWgt.getSelection()){
					weightstate = true;
					rdbtnWgt.setSelection(false);
					String msg = "Multiple weight sets take precidence over multiple weights. Multiple weights have been shut off.";
					msg += "\nUncheck Multiple weight sets to restore Multiple weights state.";
					//JOptionPane.showMessageDialog(null, msg, "Warning", JOptionPane.WARNING_MESSAGE);
				}
				else {
					weightstate = false;
				}
				*/
				
				if (rdbtnWgt.getSelection()){
					weightstate = true;
					MessageBox messageBox = new MessageBox(shShell, SWT.ICON_WARNING | SWT.OK );
		        
			        messageBox.setText("Warning");
					String msg = "Multiple weight sets take precidence over multiple weights. Multiple weights have been shut off.";
					msg += "\nUncheck Multiple weight sets to restore Multiple weights state.";
			        messageBox.setMessage(msg);
			        int buttonID = messageBox.open();
			        }
				else {
					weightstate = false;
				}

				rdbtnWgt.setSelection(true);
				txtWgtsfile.setEnabled(true);
				
			}
		}	
	}
	
	/**
	 * Launch the application.
	 */
	public static void main(final String[] args) {

		//Font btnFont = new Font("Arial", Font.BOLD, 13);
		//Color bkgndcolor = new Color(204, 255, 255);
		CliquePrototypeInterface cliquedata = new CliquePrototypeInterface();

		display = new Display();
		shShell = new Shell(display);
		shShell.setText("Clique");
        RowLayout rowLayout = new RowLayout();
        shShell.setLayout(rowLayout);

		shShell.setBounds(100, 100, 700, 350);
		//Device device = Display.getCurrent ();
		
		//shShell.setBackgroundColor(bkgndcolor);
		//shShell.setBackground(new Color(204, 255, 255));

		btnInputFile = new Button(shShell, SWT.PUSH);
		btnInputFile.setText("Input File");
	    btnInputFile.setLayoutData(new RowData(180, 28));
		btnInputFile.addSelectionListener(new SelectionAdapter() {
			public void widgetSelected(SelectionEvent event) {
				FileDialog fd = new FileDialog(shShell, SWT.OPEN);
				String selected = fd.open();
				if (selected != null)
				{
					txtInputFile.setText(selected);
				}
   			}
		});
	 
	    txtInputFile = new Text(shShell, SWT.SINGLE | SWT.BORDER);
		txtInputFile.setLayoutData(new RowData(500, 20));
		txtInputFile.setText("infile");

		
		btnOutputFile = new Button(shShell, SWT.PUSH);
	    btnOutputFile.setText("Output File");
        btnOutputFile.setLayoutData(new RowData(180, 28));
	    btnOutputFile.addSelectionListener(new SelectionAdapter() {
			public void widgetSelected(SelectionEvent event) {
				FileDialog fd = new FileDialog(shShell, SWT.OPEN);
				String selected = fd.open();
				if (selected != null)
				{
					txtOutputFile.setText(selected);
				}
	   		}
		});

        txtOutputFile = new Text(shShell, SWT.SINGLE | SWT.BORDER);
        txtOutputFile.setLayoutData(new RowData(500, 20));
	    txtOutputFile.setText("outfile");
	    
		
	    btnOutTree = new Button(shShell, SWT.PUSH);
	    btnOutTree.setText("Output Tree");
	    btnOutTree.setLayoutData(new RowData(180, 28));
	    btnOutTree.addSelectionListener(new SelectionAdapter() {
			public void widgetSelected(SelectionEvent event) {
				FileDialog fd = new FileDialog(shShell, SWT.OPEN);
				String selected = fd.open();
				if (selected != null)
				{
					txtOuttree.setText(selected);
				}
	   		}
		});

        txtOuttree = new Text(shShell, SWT.SINGLE | SWT.BORDER);
        txtOuttree.setLayoutData(new RowData(500, 20));
        txtOuttree.setText("outtree");
        
		
		rdbtnFac = new Button(shShell,SWT.PUSH);
		rdbtnFac.setText("Use Factors File");
		rdbtnFac.setLayoutData(new RowData(180, 28));
		rdbtnFac.addSelectionListener(new SelectionAdapter() {
			public void widgetSelected(SelectionEvent event) {
				FacToggle();
	   		}
		});

		txtFactfile = new Text(shShell, SWT.SINGLE | SWT.BORDER);
		txtFactfile.setLayoutData(new RowData(500, 20));
		txtFactfile.setText("factfile");
		txtFactfile.setEnabled(false);
		

		rdbtnWgt = new Button(shShell,SWT.PUSH);
		rdbtnWgt.setText("Use Site Weights File");
		rdbtnWgt.setLayoutData(new RowData(180, 26));
		//rdbtnWgt.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnWgt.addSelectionListener(new SelectionAdapter() {
			public void widgetSelected(SelectionEvent event) {
				WgtToggle();
	   		}
		});
		txtWgtsfile = new Text(shShell, SWT.SINGLE | SWT.BORDER);
		txtWgtsfile.setLayoutData(new RowData(500, 20));
		txtWgtsfile.setText("weightfile");
		txtWgtsfile.setEnabled(false);
		

		rdbtnAnc = new Button(shShell,SWT.PUSH);
		rdbtnAnc.setText("Use Ancestral States File");
		rdbtnAnc.setLayoutData(new RowData(180, 26));
		//rdbtnAnc.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnAnc.addSelectionListener(new SelectionAdapter() {
			public void widgetSelected(SelectionEvent event) {
				AncToggle();
			}
		});
		txtAncfile = new Text(shShell, SWT.SINGLE | SWT.BORDER);
		txtAncfile.setLayoutData(new RowData(500, 20));
		txtAncfile.setText("ancfile");
		txtAncfile.setEnabled(false);



		rdbtnCliq = new Button(shShell,SWT.PUSH);
		rdbtnCliq.setText("Minimum clique size");
		rdbtnCliq.setLayoutData(new RowData(180, 26));
		rdbtnCliq.setSelection(false);
		//rdbtnCliq.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnCliq.addSelectionListener(new SelectionAdapter() {
			public void widgetSelected(SelectionEvent event) {
				CliqSizeToggle();
	   		}
		});
		
	
		lblSize = new Label(shShell, SWT.NONE);
		lblSize.setText("Size:");
		lblSize.setLayoutData(new RowData(30, 20));
		//lblSize.setFont(new Font("Arial", Font.BOLD, 13));
		
		txtCliqueSize = new Text(shShell, SWT.SINGLE | SWT.BORDER);
		txtCliqueSize.setEnabled(false);
		txtCliqueSize.setText("1");
		txtCliqueSize.setLayoutData(new RowData(60, 20));

		
		
		rdbtnOgrp = new Button(shShell,SWT.PUSH);
		rdbtnOgrp.setText("Outgroup root");
		rdbtnOgrp.setLayoutData(new RowData(120, 26));
		rdbtnOgrp.setSelection(false);
		rdbtnOgrp.addSelectionListener(new SelectionAdapter() {
			public void widgetSelected(SelectionEvent event) {
				OutGroupToggle();
			}
		});

		
		rdbtnOutGroup = new Button(shShell, SWT.PUSH);
		rdbtnOutGroup.setEnabled(false);
		rdbtnOutGroup.setText("Outgroup species:");
		outspecies = false;
		rdbtnOutGroup.setLayoutData(new RowData(150, 26));
		rdbtnOutGroup.addSelectionListener(new SelectionAdapter() {
			public void widgetSelected(SelectionEvent event) {
				SpeciesToggle();
			}
		});

		txtOutRoot = new Text(shShell, SWT.SINGLE | SWT.BORDER);
		txtOutRoot.setLayoutData(new RowData(100, 20));
		txtOutRoot.setText("1");
		txtOutRoot.setEnabled(false);

		lblMultipleDatasets = new Label(shShell, SWT.NONE);
		lblMultipleDatasets.setText("Multiple sets in Input file:");
		lblMultipleDatasets.setLayoutData(new RowData(140, 20));
	
		Group WgtGrp = new Group(shShell, SWT.NONE);
		WgtGrp.setLayout(new RowLayout(SWT.HORIZONTAL));
		rdbtnMdsN = new Button(WgtGrp, SWT.RADIO);
		rdbtnMdsN.setEnabled(true);
		rdbtnMdsN.setSelection(true);
		rdbtnMdsN.setText("No");
		rdbtnMdsN.setLayoutData(new RowData(50, 26));
		rdbtnMdsN.addSelectionListener(new SelectionAdapter() {
			public void widgetSelected(SelectionEvent event) {
				MultiSetToggle(0); //Wgttype.NONE);
			}
		});

		rdbtnMdsY = new Button(WgtGrp, SWT.RADIO);
		rdbtnMdsY.setEnabled(true);
		rdbtnMdsY.setSelection(false);
		rdbtnMdsY.setText("Data sets");
		rdbtnMdsY.setLayoutData(new RowData(80, 26));
		rdbtnMdsY.addSelectionListener(new SelectionAdapter() {
			public void widgetSelected(SelectionEvent event) {
				MultiSetToggle(1); //Wgttype.DATA);
			}
		});

		rdbtnMwsY = new Button(WgtGrp, SWT.RADIO);
		rdbtnMwsY.setEnabled(true);
		rdbtnMwsY.setSelection(false);
		rdbtnMwsY.setText("Weight sets");
		rdbtnMwsY.setLayoutData(new RowData(120, 26));
		rdbtnMwsY.addSelectionListener(new SelectionAdapter() {
			public void widgetSelected(SelectionEvent event) {
				MultiSetToggle(2); //Wgttype.WEIGHT);
			}
		});
		
		lblNumSets = new Label(shShell, SWT.NONE);
		lblNumSets.setText("Number of sets:");
		lblNumSets.setLayoutData(new RowData(100, 20));
		
		txtSets = new Text(shShell, SWT.SINGLE | SWT.BORDER);
		txtSets.setLayoutData(new RowData(100, 20));
		txtSets.setText("1");
		txtSets.setEnabled(false);

		rdbtnBrun = new Button(shShell, SWT.CHECK);
		rdbtnBrun.setSelection(false);
		rdbtnBrun.setText("Print data before run");
		rdbtnBrun.setLayoutData(new RowData(160, 26));
				
		rdbtnCmat = new Button(shShell, SWT.CHECK);
		rdbtnCmat.setSelection(false);
		rdbtnCmat.setText("Print compatibility matrix");
		rdbtnCmat.setLayoutData(new RowData(200, 26));
				
		rdbtnPtree = new Button(shShell, SWT.CHECK);
		rdbtnPtree.setSelection(true);
		rdbtnPtree.setText("Print tree");
		rdbtnPtree.setLayoutData(new RowData(100, 26));
		
		rdbtnOtree = new Button(shShell, SWT.CHECK);
		rdbtnOtree.setSelection(true);
		rdbtnOtree.setText("Write tree to Output tree");
		rdbtnOtree.setLayoutData(new RowData(200, 26));

		rbtnExecute = new Button(shShell,SWT.PUSH);
		rbtnExecute.setText("Execute");
		rbtnExecute.setLayoutData(new RowData(100, 26));
		rbtnExecute.setSelection(false);
		rbtnExecute.addSelectionListener(new SelectionAdapter() {
			public void widgetSelected(SelectionEvent event) {
				
				cliquedata.infile = txtInputFile.getText();
				cliquedata.outfile = txtOutputFile.getText();
				cliquedata.outfileopt = "w";
				cliquedata.outtree = txtOuttree.getText();
				cliquedata.outtreeopt = "w";
				cliquedata.factfile = txtFactfile.getText();
				cliquedata.ancfile = txtAncfile.getText();
				cliquedata.weightfile = txtWgtsfile.getText();
				cliquedata.useAncStates = txtAncfile.getEnabled();			
				cliquedata.useFactors = txtFactfile.getEnabled();			
				cliquedata.useWeights = txtWgtsfile.getEnabled();			
				cliquedata.numWeights = Integer.parseInt(txtSets.getText());					
				cliquedata.noCliqueSize = rdbtnCliq.getEnabled();				
				cliquedata.cliqueSize = Integer.parseInt(txtCliqueSize.getText());
				cliquedata.useOutGroupRoot = rdbtnOgrp.getEnabled();			
				cliquedata.outRoot = Integer.parseInt(txtOutRoot.getText());				
				cliquedata.multipleDataSets = rdbtnMdsY.getSelection();
				cliquedata.multipleWeightSets = rdbtnMwsY.getSelection();
				cliquedata.printRawData = rdbtnBrun.getSelection();
				cliquedata.progressBar = rdbtnPrintInd.getSelection();
				cliquedata.printMatrix = rdbtnCmat.getSelection();
				cliquedata.printTree = rdbtnPtree.getSelection();
				cliquedata.outputTree = rdbtnOtree.getSelection();
				if (checkInputVals(cliquedata, shShell))
				{
			  	    Thread cliqueThread = new Thread() {
						public void run() {
			  	   		Clique clique = (Clique) Native.loadLibrary("clique", Clique.class);
			  	        clique.clique( 
			  	        		cliquedata.infile,
			  	        		cliquedata.outfile,
			  	        		cliquedata.outfileopt,
			  	        		cliquedata.outtree,
			  	        		cliquedata.outtreeopt,
			  	        		cliquedata.factfile,
			  	        		cliquedata.ancfile,
			  	          	cliquedata.weightfile,
			  	          	cliquedata.useAncStates,
			  	        		cliquedata.useFactors,
			  	        		cliquedata.useWeights,
			  	        		cliquedata.numWeights,
			  	        		cliquedata.noCliqueSize,
			  	        		cliquedata.cliqueSize,
			  	        		cliquedata.multipleDataSets,
			  	        		cliquedata.multipleWeightSets,
			  	        		cliquedata.useOutGroupRoot,
			  	        		cliquedata.outRoot,
			  	        		cliquedata.printRawData,
			  	        		cliquedata.progressBar,
			  	        		cliquedata.printMatrix,
			  	        		cliquedata.printTree,
			  	        		cliquedata.outputTree);		
						}
			  	    };
			  	    cliqueThread.start();
				}
			}
		});

		rbtnQuit = new Button(shShell,SWT.PUSH);
		rbtnQuit.setText("Quit");
		rbtnQuit.setLayoutData(new RowData(100, 26));
		rbtnQuit.setSelection(false);
		rbtnQuit.addSelectionListener(new SelectionAdapter() {
			public void widgetSelected(SelectionEvent event) {
				System.exit(0);
			}
		});
		
		rdbtnPrintInd = new Button(shShell, SWT.CHECK);
		rdbtnPrintInd.setSelection(true);
		rdbtnPrintInd.setText("Display progress");
		rdbtnPrintInd.setLayoutData(new RowData(140, 26));
		shShell.open();

		// run the event loop as long as the window is open
		while (!shShell.isDisposed()) {
		    // read the next OS event queue and transfer it to a SWT event
		    if (!display.readAndDispatch())
		     {
		    	
			    	// if there are currently no other OS event to process
			    // sleep until the next OS event is available
		        display.sleep();
		     }
		}

		// disposes all associated windows and their components
		display.dispose();
	}

	protected static boolean checkInputVals(CliquePrototypeInterface clique, Shell shShell){
		TestFileNames test = new TestFileNames();
		
		if (!test.DuplicateFileNames(clique.infile, "Input", clique.outfile, "Output", shShell))
		{			
			return false;		
		}

		if (!test.FileAvailable(clique.infile, "Input", shShell))
		{
			return false;
		}
		
		String opt = test.FileAlreadyExists(clique.outfile, "Outfile", shShell);
		if (opt == "q")
		{
			return false;
		}
		else
		{
			clique.outfileopt = opt;
		}
		
		if (clique.printTree){
			opt = test.FileAlreadyExists(clique.outtree, "Output Tree", shShell);
			if (opt == "q")
			{
				return false;
			}
			else
			{
				clique.outtreeopt = opt;
			}
		}

		if(clique.useAncStates){
			if (!test.FileAvailable(clique.ancfile, "Ancestral States", shShell))
			{
				return false;
			}
		}
		
		if(clique.useFactors){
			if (!test.FileAvailable(clique.factfile, "Factors", shShell))
			{
				return false;
			}
		}
		
		if(clique.useWeights){
			if (!test.FileAvailable(clique.weightfile, "Weights", shShell))
			{
				return false;
			}
		}
		return true;
	}
	protected void runCliqueThreads() { 
    	try
    	{
    		// see if library exists
    		Native.loadLibrary("clique", Clique.class);
		}
    	catch(UnsatisfiedLinkError e)
    	{
     		new TestFileNames().LibraryMissing("Clique", shShell);
    		return;
    	}
		try 
		{
	  	    Thread cliqueRunThread = new Thread() {
	  	    	public void run() {
		  	    	  
		  			// at this point we hook into the C code
		  	   		Clique clique = (Clique) Native.loadLibrary("clique", Clique.class);
		  	        clique.clique(
		  	        		infile,
		  	        		outfile,
		  	        		outfileopt,
		  	        		outtree,
		  	        		outtreeopt,
		  	        		factfile,
		  	        		ancfile,
	  	          		weightfile,
	  	          		useAncStates,
		  	        		useFactors,
		  	        		useWeights,
		  	        		numWeights,
		  	        		noCliqueSize,
		  	        		cliqueSize,
		  	        		multipleDataSets,
		  	        		multipleWeightSets,
		  	        		useOutGroupRoot,
		  	        		outRoot,
		  	        		printRawData,
		  	        		progressBar,
		  	        		printMatrix,
		  	        		printTree,
		  	        		outputTree);		
					}
		  	    };
	  	  	cliqueRunThread.start();
	  	  	if (progressBar)
	  	  	{
		  	  	do
		  	  	{
					Thread.sleep(1000);
					updateProgress();
				} while (cliqueRunThread.isAlive());
	  	  	}
		} 
		catch (InterruptedException e) 
		{
	  	  	if (progressBar)
	  	  	{
	  	  		updateProgress();
	  	  	}
		}
	}
	 
	private void updateProgress(){
		display.getDefault().asyncExec(new Runnable(){
			public void run() {
				MessageBox messageBox = new MessageBox(shShell, SWT.ICON_INFORMATION | SWT.OK );	        
		        messageBox.setText("Test");
		        int buttonID = messageBox.open();
		        switch (buttonID) {
		        case SWT.OK:
		        }
				/*
				DisplayProgress oldp = dp;
				dp = new DisplayProgress(inTitle,inCurdir);
				// doing the dispose in this order keeps the main window from flickering
				if (oldp != null) 
				{
					oldp.dispose();
				}
				*/
			}
		});
	}
	
}
