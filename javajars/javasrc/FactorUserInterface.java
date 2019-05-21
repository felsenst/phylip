package phylip;

import java.awt.EventQueue;
import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;

import javax.swing.JFrame;
import javax.swing.JButton;
import javax.swing.JTextField;
import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;
import javax.swing.JFileChooser;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JRadioButton;
import javax.swing.SwingConstants;
import javax.swing.SwingUtilities;

import utilities.DisplayProgress;
import utilities.TestFileNames;

import com.sun.jna.Library;
import com.sun.jna.Native;

import java.awt.Font;
import java.awt.Color;

import net.miginfocom.swing.MigLayout;
import javax.swing.JScrollPane;
import java.awt.Dimension;
import javax.swing.JPanel;

public class FactorUserInterface {
   public interface Factor extends Library {
        public void  factor(
        		String infile,
        		String outfile,
        		String ancfile,
        		String factfile,
        		boolean outputAncStates,
        		boolean outputFactors,
        		boolean printInd,
        		String outfileopt);	
    }

	public class FactorData{
		String infile;
		String outfile;
		String ancfile;
		String factfile;
		boolean outputAncStates;
		boolean outputFactors;
		boolean printInd;
		String outfileopt;
	}
	
	private FactorData inputvals;
	private String inTitle;
	private String inCurdir;
	private DisplayProgress dp;
	private String filedir;
	private boolean phylipCall;
	
	private JFrame frmFactorControls;
	private JButton btnInputFile;
	private JTextField txtInputFile;
	private JButton btnOutFile;
	private JTextField txtOutFile;
	private JButton btnAncStateFile;
	private JTextField txtAncfile;
	private JButton btnFactorFile;
	private JTextField txtFactfile;
	private JLabel lblUseAncStates;
	private JRadioButton rdbtnAncY;
	private JRadioButton rdbtnAncN;
	private JLabel lblUseFactors;
	private JRadioButton rdbtnFactYes;
	private JRadioButton rdbtnFactNo;
	private JButton btnQuit;
	private JButton btnExecute;
	private JLabel lblPrintInd;
	private JRadioButton rbtnPrintIndYes;
	private JRadioButton rbtnPrintIndNo;
	
	private JScrollPane scrollPane;
	private JPanel panel;


	/**
	 * Launch the application.
	 */
	public static void main(final String[] args) {
		EventQueue.invokeLater(new Runnable() {
			public void run() {
				try {
					FactorUserInterface window = new FactorUserInterface(args);
					window.frmFactorControls.setVisible(true);
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
		});
	}
	
	protected void ChooseFile(JTextField file) {		
		//Construct a new file choose whose default path is the path to this executable, which 
		//is returned by System.getProperty("user.dir")
		
		JFileChooser fileChooser = new JFileChooser(filedir);

		int option = fileChooser.showOpenDialog(frmFactorControls.getRootPane());
		if (option == JFileChooser.APPROVE_OPTION) {
			File selectedFile = fileChooser.getSelectedFile();
			filedir = fileChooser.getCurrentDirectory().getAbsolutePath();
			file.setText(selectedFile.getPath());
		}	
	}
	
	/**
	 * Create the application.
	 */
	public FactorUserInterface(String[] args) {
		phylipCall = false;
		if (args.length > 0)
		{
			if (args[0].contains("Phylip"))
			{
				phylipCall = true;
			}
		}
		initialize();
	}
	// event handlers

	protected void AncToggle(boolean isyes) {		
		if (isyes){
			rdbtnAncY.setSelected(true);
			rdbtnAncN.setSelected(false);
			btnAncStateFile.setEnabled(true);
			txtAncfile.setEnabled(true);
		}	
		else{
			rdbtnAncY.setSelected(false);
			rdbtnAncN.setSelected(true);
			btnAncStateFile.setEnabled(false);
			txtAncfile.setEnabled(false);
		}		
	}
	
	protected void FacToggle(boolean isyes) {		
		if (isyes){
			rdbtnFactYes.setSelected(true);
			rdbtnFactNo.setSelected(false);
			btnFactorFile.setEnabled(true);
			txtFactfile.setEnabled(true);
		}	
		else{
			rdbtnFactYes.setSelected(false);
			rdbtnFactNo.setSelected(true);
			btnFactorFile.setEnabled(false);
			txtFactfile.setEnabled(false);
		}		
	}

	protected void PrintIndToggle(boolean isPrintInd) {
		if (isPrintInd) {
			rbtnPrintIndYes.setSelected(true);
			rbtnPrintIndNo.setSelected(false);
		} else {
			rbtnPrintIndYes.setSelected(false);
			rbtnPrintIndNo.setSelected(true);
		}
	}

	/**
	 * Initialize the contents of the frame.
	 */
	private void initialize() {
		
		filedir = System.getProperty("user.dir");
		
		frmFactorControls = new JFrame();
		frmFactorControls.getContentPane().setBackground(new Color(204, 255, 255));
		frmFactorControls.setFont(new Font("Arial", Font.BOLD, 13));
		frmFactorControls.setBackground(new Color(204, 255, 255));
		frmFactorControls.setTitle("Factor");
		frmFactorControls.setBounds(100, 100, 500, 300);
		frmFactorControls.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frmFactorControls.setPreferredSize(new Dimension(frmFactorControls.getBounds().width, frmFactorControls.getBounds().height));
		
		scrollPane = new JScrollPane();
		scrollPane.setPreferredSize(frmFactorControls.getPreferredSize());
		frmFactorControls.getContentPane().add(scrollPane);
		
		panel = new JPanel();
		panel.setPreferredSize(frmFactorControls.getPreferredSize());
		scrollPane.setViewportView(panel);
		panel.setLayout(new MigLayout("", "[pref!,grow][30.00,grow][pref!,grow]", "[][][][]"));
		
		btnInputFile = new JButton("Input File");
		btnInputFile.setFont(new Font("Arial", Font.PLAIN, 13));
		btnInputFile.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				ChooseFile(txtInputFile);
			}
		});
		panel.add(btnInputFile, "cell 0 0,growx");
		
		txtInputFile = new JTextField();
		txtInputFile.setText("infile");
		txtInputFile.setFont(new Font("Arial", Font.PLAIN, 13));
		panel.add(txtInputFile, "cell 1 0 2 1,growx");
		
		btnOutFile = new JButton("Output File");
		btnOutFile.setFont(new Font("Arial", Font.PLAIN, 13));
		btnOutFile.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				ChooseFile(txtOutFile);
			}
		});
		panel.add(btnOutFile, "cell 0 1,growx");
		
		txtOutFile = new JTextField();
		txtOutFile.setText("outfile");
		txtOutFile.setFont(new Font("Arial", Font.PLAIN, 13));
		panel.add(txtOutFile, "cell 1 1 2 1,growx");
		
		btnAncStateFile = new JButton("Ancestral States File");
		btnAncStateFile.setFont(new Font("Arial", Font.PLAIN, 13));
		btnAncStateFile.setEnabled(false);
		btnAncStateFile.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				ChooseFile(txtAncfile);
			}
		});
		panel.add(btnAncStateFile, "cell 0 2,growx");
		
		txtAncfile = new JTextField();
		txtAncfile.setEnabled(false);
		txtAncfile.setText("ancfile");
		txtAncfile.setFont(new Font("Arial", Font.PLAIN, 13));
		panel.add(txtAncfile, "cell 1 2 2 1,growx");
		
		btnFactorFile = new JButton("Factors File");
		btnFactorFile.setFont(new Font("Arial", Font.BOLD, 13));
		btnFactorFile.setEnabled(false);
		btnFactorFile.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				ChooseFile(txtFactfile);
			}
		});
		panel.add(btnFactorFile, "cell 0 3,growx");
		
		txtFactfile = new JTextField();
		txtFactfile.setFont(new Font("Arial", Font.PLAIN, 13));
		txtFactfile.setEnabled(false);
		txtFactfile.setText("factfile");
		panel.add(txtFactfile, "cell 1 3 2 1,growx");
		
		lblUseAncStates = new JLabel("Output ancestral states:");
		lblUseAncStates.setHorizontalAlignment(SwingConstants.RIGHT);
		lblUseAncStates.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblUseAncStates, "flowx,cell 0 4 2 1,alignx right");
		
		rdbtnAncY = new JRadioButton("Yes");
		rdbtnAncY.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnAncY.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				AncToggle(true);				
			}
		});
		rdbtnAncY.setSelected(false);
		panel.add(rdbtnAncY, "cell 2 4");
		
		rdbtnAncN = new JRadioButton("No");
		rdbtnAncN.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnAncN.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				AncToggle(false);				
			}
		});
		rdbtnAncN.setSelected(true);
		panel.add(rdbtnAncN, "cell 2 4");
		
		
		lblUseFactors = new JLabel("Output factors:");
		lblUseFactors.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblUseFactors, "flowx,cell 0 5 2 1,alignx right");
		
		rdbtnFactYes = new JRadioButton("Yes");
		rdbtnFactYes.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnFactYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				FacToggle(true);				
			}
		});
		rdbtnFactYes.setSelected(false);
		panel.add(rdbtnFactYes, "cell 2 5");
		
		rdbtnFactNo = new JRadioButton("No");
		rdbtnFactNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnFactNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				FacToggle(false);				
			}
		});
		rdbtnFactNo.setSelected(true);
		panel.add(rdbtnFactNo, "cell 2 5");
		

		lblPrintInd = new JLabel("Display progress:");
		lblPrintInd.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblPrintInd, "flowx,cell 0 6 2 1,alignx right");

		rbtnPrintIndYes = new JRadioButton("Yes");
		rbtnPrintIndYes.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnPrintIndYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				PrintIndToggle(true);
			}
		});
		rbtnPrintIndYes.setSelected(true);
		panel.add(rbtnPrintIndYes, "cell 2 6");

		rbtnPrintIndNo = new JRadioButton("No");
		rbtnPrintIndNo.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnPrintIndNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				PrintIndToggle(false);
			}
		});
		rbtnPrintIndNo.setSelected(false);
		panel.add(rbtnPrintIndNo, "cell 2 6");
		
		btnExecute = new JButton("Execute");
		btnExecute.setFont(new Font("Arial", Font.BOLD, 13));
		btnExecute.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				inputvals = new FactorData();
				inputvals.infile = txtInputFile.getText();
				inputvals.outfile = txtOutFile.getText();
				inputvals.ancfile = txtAncfile.getText();
				inputvals.factfile = txtFactfile.getText();
				inputvals.outfileopt = "w";
				inputvals.outputAncStates = rdbtnAncY.isSelected();
				inputvals.outputFactors = rdbtnFactYes.isSelected();
				inputvals.printInd = rbtnPrintIndYes.isSelected();
					
				btnExecute.setEnabled(false);	
				String title = "Factor Progress";
				String curdir = System.getProperty("user.dir");
				curdir += "/progress.txt";
				File fl = new File(curdir);
				fl.delete();
				inTitle = title;
				inCurdir = curdir;
	
				if (checkInputVals())
				{
					
			  	    Thread factorThread = new Thread() {
						public void run() {
							runFactorThreads();
						}
			  	    };
			  	  factorThread.start();
				}
				btnExecute.setEnabled(true);
			}
		});
		panel.add(btnExecute, "cell 2 7,growx");
		
		btnQuit = new JButton("Quit");
		btnQuit.setFont(new Font("Arial", Font.BOLD, 13));
		btnQuit.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				if(phylipCall)
				{
					frmFactorControls.dispose();
					return;
				}
				else
				{
					System.exit(0);
				}
			}
		});
		panel.add(btnQuit, "cell 2 7,growx");
	}
	
	public boolean checkInputVals(){
		TestFileNames test = new TestFileNames();
		
		if (!test.DuplicateFileNames(inputvals.infile, "Input", inputvals.outfile, "Output"))
		{			
			return false;		
		}

		if (!test.FileAvailable(inputvals.infile, "Input"))
		{
			return false;
		}
		else
		{
			if (!CheckInfile(inputvals.infile))
			{
				return false;
			}
		}
		
		String opt = test.FileAlreadyExists(inputvals.outfile, "Outfile");
		if (opt == "q")
		{
			return false;
		}
		else
		{
			inputvals.outfileopt = opt;
		}
		
		if (inputvals.outputAncStates)
		{ 
			if (!test.FileAvailable(inputvals.ancfile, "Ancestral State"))
			{
				return false;
			}
		}

		if (inputvals.outputFactors)
		{
			if (!test.FileAvailable(inputvals.factfile, "Factor"))
			{
				return false;
			}
		}
		return true;
	}
		
	protected void runFactorThreads() {
	  	try
	  	{
    		// see if library exists
    		Native.loadLibrary("factor", Factor.class);
	  	}
    	catch(UnsatisfiedLinkError e)
    	{
     		new TestFileNames().LibraryMissing("Factor");
    		return;
    	}
		try 
		{
	  	    Thread dollopRunThread = new Thread() {
		  	      public void run() {
		  	    	  
		  			// at this point we hook into the C code
		  			Factor factor = (Factor) Native.loadLibrary("factor", Factor.class);
		  	        factor.factor(
		  	        		inputvals.infile,
		  	        		inputvals.outfile,
		  	        		inputvals.ancfile,
		  	           		inputvals.factfile,
		  	          		inputvals.outputAncStates,
		  	        		inputvals.outputFactors,
		  	        		inputvals.printInd,
		  	           	    inputvals.outfileopt);	
			  	    };
	  	    };
	  	    dollopRunThread.start();

	  	    if (inputvals.printInd)
	  	    {
		  	  	do
		  	  	{
					Thread.sleep(1000);
					updateProgress();
				} while (dollopRunThread.isAlive());
	  	    }
		} 
		catch (InterruptedException e) 
		{
			if (inputvals.printInd)
			{
				updateProgress();
			}
		}
	}
		  
	private void updateProgress(){
		SwingUtilities.invokeLater(new Runnable(){
			public void run() {
				DisplayProgress oldp = dp;
				dp = new DisplayProgress(inTitle,inCurdir);
				// doing the dispose in this order keeps the main window from flickering
				if (oldp != null) 
				{
					oldp.dispose();
				}
			}
		});
	}
	
	public boolean CheckInfile(String infile)
	{
		Scanner scanfile;
		try {
			scanfile = new Scanner(new File(infile));
			if (scanfile.hasNextLine())
			{
				// throw away first line, we aren't testing it
				scanfile.nextLine();
			}
			else
			{
				String msg = "Input file: ";
				msg += infile;
				msg += " is empty.";
				JOptionPane.showMessageDialog(null, msg, "Error", JOptionPane.ERROR_MESSAGE);
				return false;			
			}
			
		    int charindex = 0;
			if (scanfile.hasNextLine())
			{
				String curline = scanfile.nextLine();
			    Scanner scanline = new Scanner(curline);
			    if (scanline.hasNextInt())
			    {
				    charindex = scanline.nextInt();
				    if (charindex != 1)
				    {
						String msg = "Input file: ";
						msg += infile;
						msg += " Character definitions start with: ";
						msg += charindex;
						msg += ". They must start with 1.";
						JOptionPane.showMessageDialog(null, msg, "Error", JOptionPane.ERROR_MESSAGE);
						return false;
				    }
			    }
			}
			else
			{
				String msg = "Input file: ";
				msg += infile;
				msg += " is disordered. Character data must follow the first line.";
				JOptionPane.showMessageDialog(null, msg, "Error", JOptionPane.ERROR_MESSAGE);
				return false;			
			}
			
		    // test character order
			int lastchar = charindex;
			while (charindex != 999)
			{		    
				String curline = scanfile.nextLine();
			    Scanner scanline = new Scanner(curline);
		    	charindex = scanline.nextInt();
		    	if (charindex < lastchar)
		    	{
					String msg = "Input file: ";
					msg += infile;
					msg += " Character: ";
					msg += charindex;
					msg += " is out of order.";
					msg += " They must be in increasing order.";
					JOptionPane.showMessageDialog(null, msg, "Error", JOptionPane.ERROR_MESSAGE);
					return false;
		    	}
		    	lastchar = charindex;
			}
	
		    
		    // test species names
		    while (scanfile.hasNextLine()) 
			{		    
				String curline = scanfile.nextLine();
			    Scanner scanline = new Scanner(curline);
			    String name = scanline.next();
			    if (name.contains("(") ||
			    	name.contains(")") ||
			    	name.contains(":") ||
			    	name.contains(";") ||
			    	name.contains(",") ||
			    	name.contains("[") ||
			    	name.contains("]") )
			    	 {
						String msg = "Input file: ";
						msg += infile;
						msg += " Species: ";
						msg += name;
						msg += " may not contain any of the following characters ( ) : ; , [ ] ";
						JOptionPane.showMessageDialog(null, msg, "Error", JOptionPane.ERROR_MESSAGE);
						return false;			    	
			    	 }
			}    
			return true;
		}
		catch (FileNotFoundException e) {
			return false;
		}
	}
}


