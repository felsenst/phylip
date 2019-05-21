package phylip;
import java.awt.EventQueue;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JRadioButton;
import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;
import java.io.File;
import javax.swing.JFileChooser;
import javax.swing.JTextField;
import javax.swing.JButton;
import javax.swing.SwingConstants;
import javax.swing.SwingUtilities;

import java.awt.Font;
import java.awt.Color;
import javax.swing.DefaultComboBoxModel;
import javax.swing.JComboBox;

import com.sun.jna.Library;
import com.sun.jna.Native;

import utilities.DisplayProgress;
import utilities.TestFileNames;

import net.miginfocom.swing.MigLayout;
import javax.swing.JScrollPane;
import java.awt.Dimension;
import javax.swing.JPanel;

public class GenDistUserInterface {
	public interface GenDist extends Library {
        public void gendist(
        		String infile,
        		String outfile,
        		String outfileopt,
        		boolean AllAlleles,
        		String DistMeas,
        		boolean SquareMat,
        		boolean AnalyzeMult,
        		int NumMult,
        		boolean PrintInd);
    }

	public class GenDistData {
		String infile;
		String outfile;
		String outfileopt;
		boolean AllAlleles;
		String DistMeas;
		boolean SquareMat;
		boolean AnalyzeMult;
		int NumMult;
		boolean PrintInd;
	}
	
	private GenDistData inputvals;
	private String inTitle;
	private String inCurdir;
	private DisplayProgress dp;
	private String filedir;
	private boolean phylipCall;

	private JFrame frmGenDistControls;
	private JTextField txtInputFile;
	private JButton btnInputFile;
	private JButton btnExecute;
	private JButton btnQuit;
	private JTextField txtOutputFile;
	private JButton btnOutputFile;
	private JLabel lblAlleles;
	private JRadioButton rbtnAllelYes;
	private JRadioButton rbtnAllelNo;
	private JLabel lblDMForm;
	private JRadioButton rbtnDMSquare;
	private JRadioButton rbtnDMLowerT;
	private JLabel lblMultData;
	private JRadioButton rbtnMultDataYes;
	private JRadioButton rbtnMultDataNo;
	private JLabel lblNSets;
	private JTextField txtNSets;
	private JLabel lblPrintInd;
	private JRadioButton rbtnPrintIndYes;
	private JRadioButton rbtnPrintIndNo;
	private JLabel lblDistMeas;
	private JComboBox cmbxDistMeas;
	
	private JScrollPane scrollPane;
	private JPanel panel;

	/**
	 * Launch the application
	 */
	public static void main(final String[] args) {
		EventQueue.invokeLater(new Runnable() {
			public void run() {
				try {
					GenDistUserInterface window = new GenDistUserInterface(args);
					window.frmGenDistControls.setVisible(true);
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
		});
	}

	protected void ChooseFile(JTextField file) {
		// Construct a new file choose whose default path is the path to this
		// executable, which
		// is returned by System.getProperty("user.dir")

		JFileChooser fileChooser = new JFileChooser(filedir);

		int option = fileChooser.showOpenDialog(frmGenDistControls.getRootPane());
		if (option == JFileChooser.APPROVE_OPTION) {
			File selectedFile = fileChooser.getSelectedFile();
			filedir = fileChooser.getCurrentDirectory().getAbsolutePath();
			file.setText(selectedFile.getPath());
		}
	}
	
	protected void AllelesToggle(boolean isYes){
		if (isYes)
		{
			rbtnAllelYes.setSelected(true);
			rbtnAllelNo.setSelected(false);
		} else {
			rbtnAllelYes.setSelected(false);
			rbtnAllelNo.setSelected(true);
		}
	}

	protected void DMSquare(boolean issquare) {
		if (issquare) {
			rbtnDMSquare.setSelected(true);
			rbtnDMLowerT.setSelected(false);
		} else {
			rbtnDMSquare.setSelected(false);
			rbtnDMLowerT.setSelected(true);
		}
	}

	protected void MultToggle(boolean isMult) {
		if (isMult) {
			rbtnMultDataYes.setSelected(true);
			rbtnMultDataNo.setSelected(false);
			lblNSets.setEnabled(true);
			txtNSets.setEnabled(true);
		} else {
			rbtnMultDataYes.setSelected(false);
			rbtnMultDataNo.setSelected(true);
			lblNSets.setEnabled(false);
			txtNSets.setEnabled(false);
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
	protected void MethodToggle(int selected){
	}

	/**
	 * Create the application.
	 */
	public GenDistUserInterface(String[] args) {
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
	/**
	 * Initialize the contents of the frame.
	 */
	
	
	private void initialize() {
		
		filedir = System.getProperty("user.dir");

		frmGenDistControls = new JFrame();
		frmGenDistControls.setBackground(new Color(204, 255, 255));
		frmGenDistControls.setTitle("Gendist");
		frmGenDistControls.setBounds(100, 100, 630, 300);
		frmGenDistControls.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frmGenDistControls.setPreferredSize(new Dimension(frmGenDistControls.getBounds().width, frmGenDistControls.getBounds().height));
		
		scrollPane = new JScrollPane();
		scrollPane.setPreferredSize(frmGenDistControls.getPreferredSize());
		frmGenDistControls.getContentPane().add(scrollPane);
		
		panel = new JPanel();
		panel.setPreferredSize(frmGenDistControls.getPreferredSize());
		scrollPane.setViewportView(panel);
		panel.setLayout(new MigLayout("", "[pref!,grow][30.00,grow][pref!,grow]", "[][][][]"));

		btnInputFile = new JButton("Input File");
		btnInputFile.setFont(new Font("Arial", Font.BOLD, 13));
		btnInputFile.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				ChooseFile(txtInputFile);
			}
		});
		panel.add(btnInputFile, "cell 0 0,growx");

		txtInputFile = new JTextField();
		txtInputFile.setFont(new Font("Arial", Font.PLAIN, 13));
		txtInputFile.setText("infile");
		txtInputFile.setBounds(166, 12, 419, 20);
		panel.add(txtInputFile, "cell 1 0 2 1,growx");
		
		btnOutputFile = new JButton("Output File");
		btnOutputFile.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				ChooseFile(txtOutputFile);
			}
		});
		btnOutputFile.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(btnOutputFile, "cell 0 1,growx");
		
		txtOutputFile = new JTextField();
		txtOutputFile.setFont(new Font("Arial", Font.PLAIN, 13));
		txtOutputFile.setText("outfile");
		panel.add(txtOutputFile, "cell 1 1 2 1,growx");
		
		lblAlleles = new JLabel("All alleles present at each locus:");
		lblAlleles.setHorizontalAlignment(SwingConstants.RIGHT);
		lblAlleles.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblAlleles, "flowx,cell 0 2 2 1,alignx right");
		
		rbtnAllelYes = new JRadioButton("Yes");
		rbtnAllelYes.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnAllelYes.setBackground(new Color(204, 255, 255));
		rbtnAllelYes.setBounds(273, 72, 56, 23);
		rbtnAllelYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				AllelesToggle(true);
			}
		});
		panel.add(rbtnAllelYes, "cell 2 2");
		
		rbtnAllelNo = new JRadioButton("No, one absent at each locus");
		rbtnAllelNo.setSelected(true);
		rbtnAllelNo.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnAllelNo.setBackground(new Color(204, 255, 255));
		rbtnAllelNo.setBounds(328, 72, 220, 23);
		rbtnAllelNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				AllelesToggle(false);
			}
		});
		panel.add(rbtnAllelNo, "cell 2 2");
		
		lblDistMeas = new JLabel("Distance Measure:");
		lblDistMeas.setHorizontalAlignment(SwingConstants.RIGHT);
		lblDistMeas.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblDistMeas, "flowx,cell 0 3 2 1,alignx right");
		
		cmbxDistMeas = new JComboBox();
		cmbxDistMeas.setModel(new DefaultComboBoxModel(new String[] {"Nei genetic distance", "Cavalli-Sforza chord measure", "Reynolds genetic distance", "Delta Mu Squared microsatellite distance"}));
		cmbxDistMeas.setSelectedIndex(0);
		cmbxDistMeas.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				MethodToggle(cmbxDistMeas.getSelectedIndex());
			}
		});
		cmbxDistMeas.setFont(new Font("Arial", Font.PLAIN, 13));
		panel.add(cmbxDistMeas, "cell 2 3,growx");
		
		lblDMForm = new JLabel("Distance Matrix Form:");
		lblDMForm.setHorizontalAlignment(SwingConstants.RIGHT);
		lblDMForm.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblDMForm, "flowx,cell 0 4 2 1,alignx right");
	
		rbtnDMSquare = new JRadioButton("Square");
		rbtnDMSquare.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				DMSquare(true);
			}
		});
		rbtnDMSquare.setSelected(true);
		rbtnDMSquare.setHorizontalAlignment(SwingConstants.LEFT);
		rbtnDMSquare.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(rbtnDMSquare, "cell 2 4");
		
		rbtnDMLowerT = new JRadioButton("Lower Triangular");
		rbtnDMLowerT.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				DMSquare(false);
			}
		});
		rbtnDMLowerT.setHorizontalAlignment(SwingConstants.LEFT);
		rbtnDMLowerT.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(rbtnDMLowerT, "cell 2 4");

		lblMultData = new JLabel("Analyze multiple data sets:");
		lblMultData.setFont(new Font("Arial", Font.BOLD, 13));
		lblMultData.setHorizontalAlignment(SwingConstants.RIGHT);
		panel.add(lblMultData, "flowx,cell 0 5 2 1,alignx right");

		rbtnMultDataYes = new JRadioButton("Yes");
		rbtnMultDataYes.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnMultDataYes.setHorizontalAlignment(SwingConstants.LEFT);
		rbtnMultDataYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				MultToggle(true);
			}
		});
		rbtnMultDataYes.setSelected(false);
		panel.add(rbtnMultDataYes, "cell 2 5");

		rbtnMultDataNo = new JRadioButton("No");
		rbtnMultDataNo.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnMultDataNo.setHorizontalAlignment(SwingConstants.LEFT);
		rbtnMultDataNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				MultToggle(false);
			}
		});
		rbtnMultDataNo.setSelected(true);
		panel.add(rbtnMultDataNo, "cell 2 5");
		
		lblNSets = new JLabel("Number:\r\n");
		lblNSets.setHorizontalAlignment(SwingConstants.RIGHT);
		lblNSets.setFont(new Font("Arial", Font.BOLD, 13));
		lblNSets.setEnabled(false);
		panel.add(lblNSets, "cell 2 5");
		
		txtNSets = new JTextField();
		txtNSets.setText("1");
		txtNSets.setFont(new Font("Arial", Font.PLAIN, 13));
		txtNSets.setEnabled(false);
		txtNSets.setColumns(5);
		panel.add(txtNSets, "cell 2 5");

		lblPrintInd = new JLabel("Display progress:");
		lblPrintInd.setFont(new Font("Arial", Font.BOLD, 13));
		lblPrintInd.setHorizontalAlignment(SwingConstants.RIGHT);
		panel.add(lblPrintInd, "flowx,cell 0 6 2 1,alignx right");

		rbtnPrintIndYes = new JRadioButton("Yes");
		rbtnPrintIndYes.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnPrintIndYes.setHorizontalAlignment(SwingConstants.LEFT);
		rbtnPrintIndYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				PrintIndToggle(true);
			}
		});
		rbtnPrintIndYes.setSelected(true);
		panel.add(rbtnPrintIndYes, "cell 2 6");

		rbtnPrintIndNo = new JRadioButton("No");
		rbtnPrintIndNo.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnPrintIndNo.setHorizontalAlignment(SwingConstants.LEFT);
		rbtnPrintIndNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				PrintIndToggle(false);
			}
		});
		rbtnPrintIndNo.setSelected(false);
		panel.add(rbtnPrintIndNo, "cell 2 6");

		btnExecute = new JButton("Execute");
		btnExecute.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				inputvals = new GenDistData();
				inputvals.infile = txtInputFile.getText();
				inputvals.outfile = txtOutputFile.getText();
				inputvals.outfileopt = "w";
				inputvals.AllAlleles = rbtnAllelYes.isSelected();
				if (cmbxDistMeas.getSelectedIndex() == 0)
				{
					inputvals.DistMeas = "Nei";
				}
				else if (cmbxDistMeas.getSelectedIndex() == 1)
				{
					inputvals.DistMeas = "CS";
				}
				else if (cmbxDistMeas.getSelectedIndex() == 2)
				{
					inputvals.DistMeas = "Reynolds";
				}
				else //if (cmbxDistMeas.getSelectedIndex() == 3)
				{
					inputvals.DistMeas = "DeltaMu";				
				}
				inputvals.SquareMat = rbtnDMSquare.isSelected();
				inputvals.AnalyzeMult = rbtnMultDataYes.isSelected();
				inputvals.NumMult =  Integer.parseInt(txtNSets.getText());;
				inputvals.PrintInd = rbtnPrintIndYes.isSelected();
				
				btnExecute.setEnabled(false);
				String title = "Gendist Progress";
				String curdir = System.getProperty("user.dir");
				curdir += "/progress.txt";
				File fl = new File(curdir);
				fl.delete();
				inTitle = title;
				inCurdir = curdir;
	
				if (checkInputVals())
				{
					
			  	    Thread genDistThread = new Thread() {
						public void run() {
							runGenDistThreads();
						}
			  	    };
			  	    genDistThread.start();
				}
				btnExecute.setEnabled(true);
			}
		});
		btnExecute.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(btnExecute, "cell 2 7,alignx center");

		btnQuit = new JButton("Quit");
		btnQuit.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				if(phylipCall)
				{
					frmGenDistControls.dispose();
					return;
				}
				else
				{
					System.exit(0);
				}
			}
		});
		btnQuit.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(btnQuit, "cell 2 7");
	}
	
	public boolean checkInputVals(){
		
		// check files
		TestFileNames test = new TestFileNames();
		
		if (!test.DuplicateFileNames(inputvals.infile, "Input", inputvals.outfile, "Output"))
		{			
			return false;		
		}

		if (!test.FileAvailable(inputvals.infile, "Input"))
		{
			return false;
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
		return true;
	}
	
	protected void runGenDistThreads() {
    	try
    	{
    		// see if library exists
    		Native.loadLibrary("gendist", GenDist.class);
		}
    	catch(UnsatisfiedLinkError e)
    	{
     		new TestFileNames().LibraryMissing("GenDist");
    		return;
    	}
		try 
		{
	  	    Thread genDistRunThread = new Thread() {
		  	      public void run() {
		  	    	  
		  			// at this point we hook into the C code
		  			GenDist gendist = (GenDist) Native.loadLibrary("gendist", GenDist.class);
		  			gendist.gendist(
		  		   		inputvals.infile,
		  		   		inputvals.outfile,
		  		   		inputvals.outfileopt,
		  		   		inputvals.AllAlleles,
		  		   		inputvals.DistMeas,
		  		    	inputvals.SquareMat,
		  		   		inputvals.AnalyzeMult,
		  		   		inputvals.NumMult,
		  		   		inputvals.PrintInd);
			  	    };
	  	    };
	  	    genDistRunThread.start();

	  	    if (inputvals.PrintInd)
	  	    {
		  	  	do
		  	  	{
					Thread.sleep(1000);
					updateProgress();
				} while (genDistRunThread.isAlive());
	  	    }
		} 
		catch (InterruptedException e) 
		{
			if (inputvals.PrintInd)
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
}
