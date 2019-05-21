package phylip;
import java.awt.EventQueue;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JRadioButton;
import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Scanner;

import javax.swing.JFileChooser;
import javax.swing.JOptionPane;
import javax.swing.JTextField;
import javax.swing.JButton;
import javax.swing.SwingConstants;
import javax.swing.SwingUtilities;

import java.awt.Font;
import java.awt.Color;
import javax.swing.DefaultComboBoxModel;
import javax.swing.JComboBox;

import utilities.DisplayProgress;
import utilities.TestFileNames;

import com.sun.jna.Library;
import com.sun.jna.Native;

import net.miginfocom.swing.MigLayout;
import javax.swing.JScrollPane;
import java.awt.Dimension;
import javax.swing.JPanel;

public class ConsenseUserInterface {
	public interface Consense extends Library {
        public void consense(
        		String intree,
        		String outfile,
        		String outfileopt,
        		String outtree,
        		String outtreeopt,
        		String ConsType,
        		double Fraction,
        		boolean OutRoot,
        		int OutNum,
        		boolean Rooted,
        		boolean PrintData,
        		boolean PrintInd,
        		boolean PrintTree,
        		boolean WriteTree);
    }
	
	public class ConsenseData {
		String intree;
		String outfile;
		String outfileopt;
		String outtree;
		String outtreeopt;
		String ConsType;
		double Fraction;
		boolean OutRoot;
		int OutNum;
		boolean Rooted;
		boolean PrintData;
		boolean PrintInd;
		boolean PrintTree;
		boolean WriteTree;
	}

	private ConsenseData inputvals;
	private String inTitle;
	private String inCurdir;
	private DisplayProgress dp;
	private String filedir;
	private boolean phylipCall;
	private boolean bootstrapCall;

	private JFrame frmConsenseControls;
	private JTextField txtInputTree;
	private JTextField txtOutputFile;
	private JButton btnInputTree;
	private JButton btnOutputFile;
	private JButton btnOutputTree;
	private JTextField txtOutputTree;
	private JLabel lblConsType;
	private JComboBox cmbxConsType;
	private JLabel lblTreeRooted;
	private JRadioButton rdbtnTreeRootedYes;
	private JRadioButton rdbtnTreeRootedNo;
	private JLabel lblOutRoot;
	private JRadioButton rdbtnOutYes;
	private JRadioButton rdbtnOutNo;
	private JTextField txtOutNumber;
	private JLabel lblOutNumber;
	private JLabel lblPrintData;
	private JRadioButton rdbtnPrintDataYes;
	private JRadioButton rdbtnPrintDataNo;
	private JLabel lblPrintInd;
	private JRadioButton rdbtnPrintIndYes;
	private JRadioButton rdbtnPrintIndNo;
	private JLabel lblPrintTree;
	private JRadioButton rdbtnPrintTreeYes;
	private JRadioButton rdbtnPrintTreeNo;
	private JLabel lblWriteTree;
	private JRadioButton rdbtnWriteTreeYes;
	private JRadioButton rdbtnWriteTreeNo;
	private JButton btnExecute;
	private JButton btnQuit;
	private JLabel lblFraction;
	private JTextField txtFraction;
	
	private JScrollPane scrollPane;
	private JPanel panel;
	private JButton btnDefaults;
	private JButton btnStored;

	/**
	 * Launch the application
	 * 	public static void main(String[] args) {.
	 */
	public static void main(final String[] args) {
		EventQueue.invokeLater(new Runnable() {
			public void run() {
				try {
					ConsenseUserInterface window = new ConsenseUserInterface(args);
					window.frmConsenseControls.setVisible(true);
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

		int option = fileChooser.showOpenDialog(frmConsenseControls.getRootPane());
		if (option == JFileChooser.APPROVE_OPTION) {
			File selectedFile = fileChooser.getSelectedFile();
			filedir = fileChooser.getCurrentDirectory().getAbsolutePath();
			file.setText(selectedFile.getPath());
		}
	}
	
	protected void TreeRootedToggle(boolean isRooted) {
		if (isRooted) {
			rdbtnTreeRootedYes.setSelected(true);
			rdbtnTreeRootedNo.setSelected(false);
			lblOutRoot.setEnabled(false);
			rdbtnOutYes.setEnabled(false);
			rdbtnOutNo.setEnabled(false);
			lblOutNumber.setEnabled(false);
			txtOutNumber.setEnabled(false);
		} else {
			rdbtnTreeRootedYes.setSelected(false);
			rdbtnTreeRootedNo.setSelected(true);
			lblOutRoot.setEnabled(true);
			rdbtnOutYes.setEnabled(true);
			rdbtnOutNo.setEnabled(true);
			OutToggle(rdbtnOutYes.isSelected());
		}
	}

	protected void OutToggle(boolean isOut) {
		if (isOut) {
			rdbtnOutYes.setSelected(true);
			rdbtnOutNo.setSelected(false);
			lblOutNumber.setEnabled(true);
			txtOutNumber.setEnabled(true);

		} else {
			rdbtnOutYes.setSelected(false);
			rdbtnOutNo.setSelected(true);
			lblOutNumber.setEnabled(false);
			txtOutNumber.setEnabled(false);
		}
	}

	protected void PrintDataToggle(boolean isPrintData) {
		if (isPrintData) {
			rdbtnPrintDataYes.setSelected(true);
			rdbtnPrintDataNo.setSelected(false);
		} else {
			rdbtnPrintDataYes.setSelected(false);
			rdbtnPrintDataNo.setSelected(true);
		}
	}

	protected void PrintIndToggle(boolean isPrintInd) {
		if (isPrintInd) {
			rdbtnPrintIndYes.setSelected(true);
			rdbtnPrintIndNo.setSelected(false);
		} else {
			rdbtnPrintIndYes.setSelected(false);
			rdbtnPrintIndNo.setSelected(true);
		}
	}

	protected void PrintTreeToggle(boolean isTree) {
		if (isTree) {
			rdbtnPrintTreeYes.setSelected(true);
			rdbtnPrintTreeNo.setSelected(false);
		} else {
			rdbtnPrintTreeYes.setSelected(false);
			rdbtnPrintTreeNo.setSelected(true);
		}
	}

	protected void WriteTreeToggle(boolean isWrite) {
		if (isWrite) {
			rdbtnWriteTreeYes.setSelected(true);
			rdbtnWriteTreeNo.setSelected(false);
			btnOutputTree.setEnabled(true);
			txtOutputTree.setEnabled(true);
		} else {
			rdbtnWriteTreeYes.setSelected(false);
			rdbtnWriteTreeNo.setSelected(true);
			btnOutputTree.setEnabled(false);
			txtOutputTree.setEnabled(false);
		}
	}
	
	protected void ConsToggle(int selected) {
		if (selected == 0)  // Majority rule (extended)
		{
			lblFraction.setEnabled(false);
			txtFraction.setEnabled(false);
		}
		else if (selected == 1) // Strict
		{
			lblFraction.setEnabled(false);
			txtFraction.setEnabled(false);
		}
		else if (selected == 2) // Majority rule
		{
			lblFraction.setEnabled(false);
			txtFraction.setEnabled(false);
		}
		else //if (selected == 3) Ml
		{
			lblFraction.setEnabled(true);
			txtFraction.setEnabled(true);
		}
	}
	
	protected void DoInit()
	{
		// reset everything if there is an init file
		getStoredSettings();
	}

	/**
	 * Create the application.
	 */
	public ConsenseUserInterface(String[] args) {
		phylipCall = false;
		bootstrapCall = false;
		if (args.length > 0)
		{
			if (args[0].contains("Phylip"))
			{
				phylipCall = true;
			}
			if (args.length > 1)
			{
				if (args[1].contains("bootstrap"))
				{
					bootstrapCall = true;
				}
			}
		}
		initialize();
		DoInit();
		
		if (bootstrapCall){
  			getDnaMLSettings();
 		}
	}
	/**
	 * Initialize the contents of the frame.
	 */
	
	
	private void initialize() {
		
		filedir = System.getProperty("user.dir");

		frmConsenseControls = new JFrame();
		frmConsenseControls.getContentPane().setBackground(new Color(204, 255, 255));
		frmConsenseControls.setBackground(new Color(204, 255, 255));
		frmConsenseControls.setTitle("Consense");
		frmConsenseControls.setFont(new Font("Arial", Font.BOLD, 13));
		frmConsenseControls.setBounds(100, 100, 600, 450);
		frmConsenseControls.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frmConsenseControls.setPreferredSize(new Dimension(frmConsenseControls.getBounds().width, frmConsenseControls.getBounds().height));
		
		scrollPane = new JScrollPane();
		scrollPane.setPreferredSize(frmConsenseControls.getPreferredSize());
		frmConsenseControls.getContentPane().add(scrollPane);
		
		panel = new JPanel();
		panel.setPreferredSize(frmConsenseControls.getPreferredSize());
		scrollPane.setViewportView(panel);
		panel.setLayout(new MigLayout("", "[pref!,grow][150px][pref!,grow]", "[][][][][][][][][][][][][]"));
		
		btnInputTree = new JButton("Input Tree");
		btnInputTree.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				ChooseFile(txtInputTree);
			}
		});
		btnInputTree.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(btnInputTree, "cell 0 0,growx");
		
		txtInputTree = new JTextField();
		txtInputTree.setFont(new Font("Arial", Font.PLAIN, 13));
		panel.add(txtInputTree, "cell 1 0 2 1,growx");
		
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
		txtOutputFile.setBounds(166, 40, 419, 20);
		panel.add(txtOutputFile, "cell 1 1 2 1,growx");
		
		btnOutputTree = new JButton("Output Tree");
		btnOutputTree.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				ChooseFile(txtOutputTree);
			}
		});
		btnOutputTree.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(btnOutputTree, "cell 0 2,growx");
		
		txtOutputTree = new JTextField();
		txtOutputTree.setFont(new Font("Arial", Font.PLAIN, 13));
		panel.add(txtOutputTree, "cell 1 2 2 1,growx");
		
		lblConsType = new JLabel("Consensus type:");
		lblConsType.setFont(new Font("Arial", Font.BOLD, 13));
		lblConsType.setHorizontalAlignment(SwingConstants.RIGHT);
		panel.add(lblConsType, "flowx,cell 0 3 2 1,alignx right");
		
		cmbxConsType = new JComboBox();
		cmbxConsType.setModel(new DefaultComboBoxModel(new String[] {"Majority rule (extended)", "Strict", "Majority rule", "Ml"}));
		cmbxConsType.setBounds(305, 112, 264, 23);
		cmbxConsType.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				ConsToggle(cmbxConsType.getSelectedIndex());
			}
		});
		panel.add(cmbxConsType, "cell 2 3,growx");

		lblFraction = new JLabel("Fraction (l) of times a branch must appear:");
		lblFraction.setHorizontalAlignment(SwingConstants.RIGHT);
		lblFraction.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblFraction, "flowx,cell 0 4 2 1,alignx right");
		
		txtFraction = new JTextField();
		txtFraction.setColumns(6);
		txtFraction.setFont(new Font("Arial", Font.PLAIN, 13));
		panel.add(txtFraction, "cell 2 4");

		lblOutRoot = new JLabel("Outgroup Root:");
		lblOutRoot.setFont(new Font("Arial", Font.BOLD, 13));
		lblOutRoot.setHorizontalAlignment(SwingConstants.RIGHT);
		panel.add(lblOutRoot, "flowx,cell 0 5 2 1,alignx right");

		rdbtnOutYes = new JRadioButton("Yes");
		rdbtnOutYes.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnOutYes.setHorizontalAlignment(SwingConstants.LEFT);
		rdbtnOutYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				OutToggle(true);
			}
		});
		panel.add(rdbtnOutYes, "cell 2 5");

		rdbtnOutNo = new JRadioButton("No, use as outgroup species");
		rdbtnOutNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnOutNo.setHorizontalAlignment(SwingConstants.LEFT);
		rdbtnOutNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				OutToggle(false);
			}
		});
		panel.add(rdbtnOutNo, "cell 2 5,growx");

		lblOutNumber = new JLabel("Number of the outgroup:");
		lblOutNumber.setFont(new Font("Arial", Font.BOLD, 13));
		lblOutNumber.setHorizontalAlignment(SwingConstants.RIGHT);
		panel.add(lblOutNumber, "flowx,cell 2 6 2 1");

		txtOutNumber = new JTextField();
		txtOutNumber.setColumns(5);
		txtOutNumber.setFont(new Font("Arial", Font.PLAIN, 13));
		panel.add(txtOutNumber, "cell 2 6");

		lblTreeRooted = new JLabel("Trees rooted:");
		lblTreeRooted.setFont(new Font("Arial", Font.BOLD, 13));
		lblTreeRooted.setHorizontalAlignment(SwingConstants.RIGHT);
		lblTreeRooted.setBounds(51, 232, 243, 14);
		panel.add(lblTreeRooted, "flowx,cell 0 7 2 1,alignx right");

		rdbtnTreeRootedYes = new JRadioButton("Yes");
		rdbtnTreeRootedYes.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnTreeRootedYes.setHorizontalAlignment(SwingConstants.LEFT);
		rdbtnTreeRootedYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				TreeRootedToggle(true);
			}
		});
		panel.add(rdbtnTreeRootedYes, "cell 2 7");

		rdbtnTreeRootedNo = new JRadioButton("No");
		rdbtnTreeRootedNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnTreeRootedNo.setHorizontalAlignment(SwingConstants.LEFT);
		rdbtnTreeRootedNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				TreeRootedToggle(false);
			}
		});
		panel.add(rdbtnTreeRootedNo, "cell 2 7");

		lblPrintData = new JLabel("Print out the sets of species:");
		lblPrintData.setFont(new Font("Arial", Font.BOLD, 13));
		lblPrintData.setHorizontalAlignment(SwingConstants.RIGHT);
		panel.add(lblPrintData, "flowx,cell 0 8 2 1,alignx right");

		rdbtnPrintDataYes = new JRadioButton("Yes");
		rdbtnPrintDataYes.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnPrintDataYes.setHorizontalAlignment(SwingConstants.LEFT);
		rdbtnPrintDataYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				PrintDataToggle(true);
			}
		});
		panel.add(rdbtnPrintDataYes, "cell 2 8");

		rdbtnPrintDataNo = new JRadioButton("No");
		rdbtnPrintDataNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnPrintDataNo.setHorizontalAlignment(SwingConstants.LEFT);
		rdbtnPrintDataNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				PrintDataToggle(false);
			}
		});
		panel.add(rdbtnPrintDataNo, "cell 2 8");

		lblPrintTree = new JLabel("Print out tree:");
		lblPrintTree.setFont(new Font("Arial", Font.BOLD, 13));
		lblPrintTree.setHorizontalAlignment(SwingConstants.RIGHT);
		panel.add(lblPrintTree, "flowx,cell 0 9 2 1,alignx right");

		rdbtnPrintTreeYes = new JRadioButton("Yes");
		rdbtnPrintTreeYes.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnPrintTreeYes.setHorizontalAlignment(SwingConstants.LEFT);
		rdbtnPrintTreeYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				PrintTreeToggle(true);
			}
		});
		panel.add(rdbtnPrintTreeYes, "cell 2 9");

		rdbtnPrintTreeNo = new JRadioButton("No");
		rdbtnPrintTreeNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnPrintTreeNo.setHorizontalAlignment(SwingConstants.LEFT);
		rdbtnPrintTreeNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				PrintTreeToggle(false);
			}
		});
		panel.add(rdbtnPrintTreeNo, "cell 2 9");

		lblWriteTree = new JLabel("Write out trees onto tree file:");
		lblWriteTree.setFont(new Font("Arial", Font.BOLD, 13));
		lblWriteTree.setHorizontalAlignment(SwingConstants.RIGHT);
		panel.add(lblWriteTree, "flowx,cell 0 10 2 1,alignx right");

		rdbtnWriteTreeYes = new JRadioButton("Yes");
		rdbtnWriteTreeYes.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnWriteTreeYes.setHorizontalAlignment(SwingConstants.LEFT);
		rdbtnWriteTreeYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				WriteTreeToggle(true);
			}
		});
		panel.add(rdbtnWriteTreeYes, "cell 2 10");

		rdbtnWriteTreeNo = new JRadioButton("No");
		rdbtnWriteTreeNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnWriteTreeNo.setHorizontalAlignment(SwingConstants.LEFT);
		rdbtnWriteTreeNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				WriteTreeToggle(false);
			}
		});
		panel.add(rdbtnWriteTreeNo, "cell 2 10");

		lblPrintInd = new JLabel("Display progress:");
		lblPrintInd.setFont(new Font("Arial", Font.BOLD, 13));
		lblPrintInd.setHorizontalAlignment(SwingConstants.RIGHT);
		panel.add(lblPrintInd, "flowx,cell 0 11 2 1,alignx right");

		rdbtnPrintIndYes = new JRadioButton("Yes");
		rdbtnPrintIndYes.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnPrintIndYes.setHorizontalAlignment(SwingConstants.LEFT);
		rdbtnPrintIndYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				PrintIndToggle(true);
			}
		});
		panel.add(rdbtnPrintIndYes, "cell 2 11");

		rdbtnPrintIndNo = new JRadioButton("No");
		rdbtnPrintIndNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnPrintIndNo.setHorizontalAlignment(SwingConstants.LEFT);
		rdbtnPrintIndNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				PrintIndToggle(false);
			}
		});
		panel.add(rdbtnPrintIndNo, "cell 2 11");
		
		btnDefaults = new JButton("Reset to Defaults");
		btnDefaults.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				resetDefaults();
			}
		});
		btnDefaults.setFont(new Font("Arial", Font.BOLD, 13));	
		panel.add(btnDefaults, "cell 0 12");
		
		btnStored = new JButton("Read Init file");
		btnStored.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				getStoredSettings();
			}
		});
		btnStored.setFont(new Font("Arial", Font.BOLD, 13));	
		panel.add(btnStored, "cell 1 12");

		btnExecute = new JButton("Execute");
		btnExecute.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				inputvals = getInputVals();
				
				btnExecute.setEnabled(false);	
				String title = "Consense Progress";
				String curdir = System.getProperty("user.dir");
				curdir += "/progress.txt";
				File fl = new File(curdir);
				fl.delete();
				inTitle = title;
				inCurdir = curdir;
	
				if (checkInputVals())
				{
					
			  	    Thread consenseThread = new Thread() {
						public void run() {
							runConsenseThreads();
						}
			  	    };
			  	    consenseThread.start();
				}
				btnExecute.setEnabled(true);
		}
		});
		btnExecute.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(btnExecute, "cell 2 12");

		btnQuit = new JButton("Quit");
		btnQuit.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				saveSettings();
				if(phylipCall)
				{
					frmConsenseControls.dispose();
					return;
				}
				else
				{
					System.exit(0);
				}
			}
		});
		btnQuit.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(btnQuit, "cell 2 12");

	}
	public boolean checkInputVals(){
		TestFileNames test = new TestFileNames();

		if(inputvals.WriteTree)
		{
			if (!test.DuplicateFileNames(inputvals.intree, "Input Tree", inputvals.outtree, "Output Tree"))
			{			
				return false;		
			}
		}

		if (!test.FileAvailable(inputvals.intree, "Input Tree"))
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
		
		if(inputvals.WriteTree)
		{
			opt = test.FileAlreadyExists(inputvals.outtree, "Output Tree");
			if (opt == "q")
			{
				return false;
			}
			else
			{
				inputvals.outtreeopt = opt;
			}
		}

		String msg;		
		if ((inputvals.Fraction < 0.5) || (inputvals.Fraction > 1.0))
		{
			msg = "Branch fraction must be between 0.5 and 1.0. ";
			JOptionPane.showMessageDialog(null, msg, "Error", JOptionPane.ERROR_MESSAGE);
			return false;
		}
		
		if (inputvals.OutNum < 1) 
		{
			msg = "Bad outgroup number. Must be greater than zero.";
			JOptionPane.showMessageDialog(null, msg, "Error", JOptionPane.ERROR_MESSAGE);
			return false;
		}
		
		return true;
	}
	
	protected void runConsenseThreads() {
    	try
    	{
    		// see if library exists
    		Native.loadLibrary("consense", Consense.class);
		}
    	catch(UnsatisfiedLinkError e)
    	{
     		new TestFileNames().LibraryMissing("Consense");
    		return;
    	}
		try 
		{
	  	    Thread consenseRunThread = new Thread() {
		  	      public void run() {
		  	    	  
		  	    	  // at this point we hook into the C code
		  	    	  Consense consense = (Consense) Native.loadLibrary("consense", Consense.class);
				      consense.consense(
				    		inputvals.intree,
				    		inputvals.outfile,
				    		inputvals.outfileopt,
				    		inputvals.outtree,
				    		inputvals.outtreeopt,
				    		inputvals.ConsType,
				    		inputvals.Fraction,
				     		inputvals.OutRoot,
				    		inputvals.OutNum,
				    		inputvals.Rooted,
				    		inputvals.PrintData,
				    		inputvals.PrintInd,
				    		inputvals.PrintTree,
				    		inputvals.WriteTree);
				    		
		  	    };
	  	    };
	  	    consenseRunThread.start();

	  	    if (inputvals.PrintInd)
	  	    {
		  	  	do
		  	  	{
					Thread.sleep(1000);
					updateProgress();
				} while (consenseRunThread.isAlive());
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
	
	protected ConsenseData getInputVals()
	{
		ConsenseData inputvals = new ConsenseData();
		
		inputvals.intree = txtInputTree.getText();
		inputvals.outfile = txtOutputFile.getText();
		inputvals.outfileopt = "w";
		inputvals.outtree = txtOutputTree.getText();
		inputvals.outtreeopt = "w";
		if (cmbxConsType.getSelectedIndex() == 0)
		{
			inputvals.ConsType = "extended";
		}
		else if (cmbxConsType.getSelectedIndex() == 1)
		{
			inputvals.ConsType = "strict";
		}
		else if (cmbxConsType.getSelectedIndex() == 2)
		{
			inputvals.ConsType = "majority";
		}
		else //cmbxConsType.getSelectedIndex() == 3)
		{
			inputvals.ConsType = "msubl";
		}
		inputvals.Fraction = Double.parseDouble(txtFraction.getText());
		inputvals.OutRoot = rdbtnOutYes.isSelected();
		inputvals.OutNum = Integer.parseInt(txtOutNumber.getText());
		inputvals.Rooted = rdbtnTreeRootedYes.isSelected();;
		inputvals.PrintData = rdbtnPrintDataYes.isSelected();
		inputvals.PrintInd = rdbtnPrintIndYes.isSelected();
		inputvals.PrintTree = rdbtnPrintTreeYes.isSelected();
		inputvals.WriteTree = rdbtnWriteTreeYes.isSelected();
		
		return inputvals;
	}
	protected void saveSettings(){
		inputvals = getInputVals();
		// there must be a better way to format this output, but this works for the prototype JRM
        try {
            BufferedWriter output = new BufferedWriter(new FileWriter("consenseInit.txt"));
    		output.write("intree : "+inputvals.intree+"\n");
    		output.write("outfile : "+inputvals.outfile+"\n");
    		//output.write("outfileopt : "+inputvals.outfileopt+"\n"); //makes no sense to save, can change between runs JRM
    		output.write("outtree : "+inputvals.outtree+"\n");
    		//output.write("outtreeopt : "+inputvals.outtreeopt+"\n"); //makes no sense to save, can change between runs JRM
    		output.write("ConsType : "+inputvals.ConsType+"\n");
    		output.write("Fraction : "+inputvals.Fraction+"\n");
     		output.write("OutRoot : "+String.format("%b",inputvals.OutRoot+"\n"));
    		output.write("OutNum : "+inputvals.OutNum+"\n");
    		output.write("Rooted : "+String.format("%b",inputvals.Rooted+"\n"));
    		output.write("PrintData : "+String.format("%b",inputvals.PrintData+"\n"));
    		output.write("PrintInd : "+String.format("%b",inputvals.PrintInd+"\n"));
    		output.write("PrintTree : "+String.format("%b",inputvals.PrintTree+"\n"));
    		output.write("WriteTree : "+String.format("%b",inputvals.WriteTree+"\n"));
        } catch ( IOException ioerr ) {
            ioerr.printStackTrace();
       }        
	}	
	
	protected void getStoredSettings(){
		// because we are setting screen values directly, this is a tedious mess to set up JRM
	    try 
	    {
	    	Scanner scanner =  new Scanner(new File("seqbootInit.txt"));
	        while (scanner.hasNextLine()){
	        	Scanner linescan =  new Scanner( scanner.nextLine());
	        	linescan.useDelimiter(" : ");
	        	String label = linescan.next();
	        	String value = linescan.next();
	       		if("intree".equals(label)){
	       			txtInputTree.setText(value);
	       		}
	       		else if("outfile".equals(label)){
	       			txtOutputFile.setText(value);
	       		}
	    		else if("outtree".equals(label)){
	    			txtOutputTree.setText(value);
	       		}
	    		else if("ConsType".equals(label)){
					if (value.contains("extended"))
					{
						cmbxConsType.setSelectedIndex(0);
					}
					else if (value.contains("strict"))
					{
						cmbxConsType.setSelectedIndex(1);				
					}
					else if (value.contains("majority"))
					{
						cmbxConsType.setSelectedIndex(2);				
					}
					else //if (value.contains("msubl"))
					{
						cmbxConsType.setSelectedIndex(3);				
					}
	       		}
	    		else if("Fraction".equals(label)){
	    			txtFraction.setText(value);
	       		}
	     		else if("OutRoot".equals(label)){
	    			if ("true".equals(value))
	    			{
	    				rdbtnOutYes.setSelected(true);
	    				rdbtnOutNo.setSelected(false);
	    			}
	    			else
	    			{
	    				rdbtnOutYes.setSelected(false);
	    				rdbtnOutNo.setSelected(true);
	    			}
	       		}
	    		else if("OutNum".equals(label)){
	    			txtOutNumber.setText(value);
	       		}
	    		else if("Rooted".equals(label)){
	    			if ("true".equals(value))
	    			{
	    				rdbtnTreeRootedYes.setSelected(true);
	    				rdbtnTreeRootedNo.setSelected(false);
	    			}
	    			else
	    			{
	    				rdbtnTreeRootedYes.setSelected(false);
	    				rdbtnTreeRootedNo.setSelected(true);
	    			}    			
	       		}
	    		else if("PrintData".equals(label)){
	    			if ("true".equals(value))
	    			{
	    				rdbtnPrintDataYes.setSelected(true);
	    				rdbtnPrintDataNo.setSelected(false);
	    			}
	    			else
	    			{
	    				rdbtnPrintDataYes.setSelected(false);
	    				rdbtnPrintDataNo.setSelected(true);
	    			}    			
	       		}
	    		else if("PrintInd".equals(label)){
	    			if ("true".equals(value))
	    			{
	    				rdbtnPrintIndYes.setSelected(true);
	    				rdbtnPrintIndNo.setSelected(false);
	    			}
	    			else
	    			{
	    				rdbtnPrintIndYes.setSelected(false);
	    				rdbtnPrintIndNo.setSelected(true);
	    			}    				       			
	       		}
	    		else if("PrintTree".equals(label)){
	    			if ("true".equals(value))
	    			{
	    				rdbtnPrintTreeYes.setSelected(true);
	    				rdbtnPrintTreeNo.setSelected(false);
	    			}
	    			else
	    			{
	    				rdbtnPrintTreeYes.setSelected(false);
	    				rdbtnPrintTreeNo.setSelected(true);
	    			}    			
	       		}
	    		else if("WriteTree".equals(label)){
	    			if ("true".equals(value))
	    			{
	    				rdbtnWriteTreeYes.setSelected(true);
	    				rdbtnWriteTreeNo.setSelected(false);
	    			}
	    			else
	    			{
	    				rdbtnWriteTreeYes.setSelected(false);
	    				rdbtnWriteTreeNo.setSelected(true);
	    			}    			
	       		}
	    		else {
	    			String msg = "Unknown label: ";
	    			msg += label;
	    			msg += " with value: ";
	    			msg += value;
	    			msg += " found in seqbootInit.txt.";
	    			JOptionPane.showMessageDialog(null, msg, "Warning", JOptionPane.WARNING_MESSAGE);			
	    		}
	        }
	    }
		catch (FileNotFoundException e)
		{
			// if it's not there, use the defaults
			resetDefaults();
		} 	
	}
	
	protected void resetDefaults()
	{
		// reset Consense defaults
		txtInputTree.setText("intree");
		txtOutputFile.setText("outfile");
		txtOutputTree.setText("outtree");
		cmbxConsType.setSelectedIndex(0);
		lblFraction.setEnabled(false);
		txtFraction.setEnabled(false);
		txtFraction.setText("0.5");
		rdbtnOutYes.setSelected(false);
		rdbtnOutNo.setSelected(true);
		lblOutNumber.setEnabled(false);
		txtOutNumber.setEnabled(false);
		txtOutNumber.setText("1");
		rdbtnTreeRootedYes.setSelected(false);
		rdbtnTreeRootedNo.setSelected(true);
		rdbtnPrintDataYes.setSelected(true);
		rdbtnPrintDataNo.setSelected(false);
		rdbtnPrintTreeYes.setSelected(true);
		rdbtnPrintTreeNo.setSelected(false);
		rdbtnWriteTreeYes.setSelected(true);
		rdbtnWriteTreeNo.setSelected(false);
		rdbtnPrintIndYes.setSelected(true);
		rdbtnPrintIndNo.setSelected(false);
	}
	
	
	protected void getDnaMLSettings()
	{
		// get the parameters needed
	    try 
	    {
	    	Scanner scanner =  new Scanner(new File("dnamlInit.txt"));
	        while (scanner.hasNextLine()){
	        	Scanner linescan =  new Scanner( scanner.nextLine());
	        	linescan.useDelimiter(" : ");
	        	String label = linescan.next();
	        	String value = linescan.next();
	      		if("outtree".equals(label)){
	      			txtInputTree.setText(value);
	        	}
	        }
	    }
		catch (FileNotFoundException e)
		{
			String msg = "Input file: dnamlInit.txt does not exist.";
			JOptionPane.showMessageDialog(null, msg, "Error", JOptionPane.ERROR_MESSAGE);			
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
