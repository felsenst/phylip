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

public class TreeDistUserInterface {
   public interface TreeDist extends Library {
        public void treedist(
        		String intree,
        		String intree2,
        		String outfile,
        		String outfileopt,
        		String DistType,
        		String DistKind,
        		boolean Rooted,
        		String OutputKind,
        		boolean PrintInd);
    }

	public class TreeDistData {
		String intree;
		String intree2;
		String outfile;
		String outfileopt;
		String DistType;
		String DistKind;
		boolean Rooted;
		String OutputKind;
		boolean PrintInd;
	}

	private TreeDistData inputvals;
	private String inTitle;
	private String inCurdir;
	private DisplayProgress dp;
	private String filedir;
	private boolean phylipCall;
	
	private String [] OutputForm1;
	private String [] OutputForm2;

	private JFrame frmTreeDistControls;
	private JTextField txtInputTree;
	private JTextField txtOutputFile;
	private JButton btnInputTree;
	private JButton btnOutputFile;
	private JLabel lblDistType;
	private JComboBox cmbxDistType;
	private JLabel lblTreeRooted;
	private JRadioButton rdbtnTreeRootedYes;
	private JRadioButton rdbtnTreeRootedNo;
	private JLabel lblPrintInd;
	private JRadioButton rdbtnPrintIndYes;
	private JRadioButton rdbtnPrintIndNo;
	private JButton btnExecute;
	private JButton btnQuit;
	private JTextField txtInputTree2;
	private JLabel lblDistKind;
	private JComboBox cmbxDistKind;
	private JButton btnInputTree2;
	private JLabel lblOutputKind;
	private JComboBox cmbxOutputKind;
	
	private JScrollPane scrollPane;
	private JPanel panel;

	/**
	 * Launch the application
	 */
	public static void main(final String[] args) {
		EventQueue.invokeLater(new Runnable() {
			public void run() {
				try {
					TreeDistUserInterface window = new TreeDistUserInterface(args);
					window.frmTreeDistControls.setVisible(true);
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

		int option = fileChooser.showOpenDialog(frmTreeDistControls.getRootPane());
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
		} else {
			rdbtnTreeRootedYes.setSelected(false);
			rdbtnTreeRootedNo.setSelected(true);
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
	
	protected void DistTypeToggle(int selected) {
		if (selected == 0)  // Branch Score
		{
		}
		else if (selected == 1) // Symmetric
		{
		}
		else //if (selected == 2) Robinson-Foulds
		{
		}
	}
	
	protected void DistKindToggle(int selected) {
		if (selected == 0)  // Adjacent pairs in single tree file, 
		{
			cmbxOutputKind.setModel(new DefaultComboBoxModel(OutputForm2));
			btnInputTree2.setEnabled(false);
			txtInputTree2.setEnabled(false);
		}
		else if (selected == 1) // All possible pairs in a single tree file
		{
			cmbxOutputKind.setModel(new DefaultComboBoxModel(OutputForm1));
			btnInputTree2.setEnabled(false);
			txtInputTree2.setEnabled(false);
		}
		else if (selected == 2) // Corresponding pairs between two tree files
		{
			cmbxOutputKind.setModel(new DefaultComboBoxModel(OutputForm2));
			btnInputTree2.setEnabled(true);
			txtInputTree2.setEnabled(true);
		}
		else //if (selected == 3) All possible pairs between two tree files
		{
			cmbxOutputKind.setModel(new DefaultComboBoxModel(OutputForm1));
			btnInputTree2.setEnabled(true);
			txtInputTree2.setEnabled(true);
		}
	}

	/**
	 * Create the application.
	 */
	public TreeDistUserInterface(String[] args) {
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
		OutputForm1 = new String[] {"Full matrix", "One pair per line, verbose", "One pair per line, sparse"};
		OutputForm2 = new String[] {"One pair per line, verbose", "One pair per line, sparse"};

		frmTreeDistControls = new JFrame();
		frmTreeDistControls.setBackground(new Color(204, 255, 255));
		frmTreeDistControls.setTitle("Treedist");
		frmTreeDistControls.setBounds(100, 100, 630, 350);
		frmTreeDistControls.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frmTreeDistControls.setPreferredSize(new Dimension(frmTreeDistControls.getBounds().width, frmTreeDistControls.getBounds().height));
		
		scrollPane = new JScrollPane();
		scrollPane.setPreferredSize(frmTreeDistControls.getPreferredSize());
		frmTreeDistControls.getContentPane().add(scrollPane);
		
		panel = new JPanel();
		panel.setPreferredSize(frmTreeDistControls.getPreferredSize());
		scrollPane.setViewportView(panel);
		panel.setLayout(new MigLayout("", "[pref!,grow][30.00,grow][pref!,grow]", "[][][][]"));
		
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
		txtInputTree.setText("intree");
		panel.add(txtInputTree, "cell 1 0 2 1,growx");
		
		btnInputTree2 = new JButton("Input Tree 2");
		btnInputTree2.setEnabled(false);
		btnInputTree2.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(btnInputTree2, "cell 0 1,growx");
		
		txtInputTree2 = new JTextField();
		txtInputTree2.setEnabled(false);
		txtInputTree2.setText("intree2");
		txtInputTree2.setFont(new Font("Arial", Font.PLAIN, 13));
		panel.add(txtInputTree2, "cell 1 1 2 1,growx");
		
		btnOutputFile = new JButton("Output File");
		btnOutputFile.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				ChooseFile(txtOutputFile);
			}
		});
		btnOutputFile.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(btnOutputFile, "cell 0 2,growx");
		
		txtOutputFile = new JTextField();
		txtOutputFile.setFont(new Font("Arial", Font.PLAIN, 13));
		txtOutputFile.setText("outfile");
		panel.add(txtOutputFile, "cell 1 2 2 1,growx");

		lblDistType = new JLabel("Distance type:");
		lblDistType.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblDistType, "flowx,cell 0 3 2 1,alignx right");
		
		cmbxDistType = new JComboBox();
		cmbxDistType.setModel(new DefaultComboBoxModel(new String[] {"Branch Score", "Symmetric", "Robinson-Foulds"}));
		cmbxDistType.setSelectedIndex(0);
		cmbxDistType.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				DistTypeToggle(cmbxDistType.getSelectedIndex());
			}
		});
		panel.add(cmbxDistType, "cell 2 3, growx");
		
		lblDistKind = new JLabel("Distances between:");
		lblDistKind.setHorizontalAlignment(SwingConstants.RIGHT);
		lblDistKind.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblDistKind, "flowx,cell 0 4 2 1,alignx right");
		
		cmbxDistKind = new JComboBox();
		cmbxDistKind.setModel(new DefaultComboBoxModel(new String[] {"Adjacent pairs in single tree file", 
																	 "All possible pairs in a single tree file", 
																	 "Corresponding pairs between two tree files", 
																	 "All possible pairs between two tree files"}));
		cmbxDistKind.setSelectedIndex(0);
		cmbxDistKind.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				DistKindToggle(cmbxDistKind.getSelectedIndex());
			}
		});
		panel.add(cmbxDistKind, "cell 2 4, growx");

		lblTreeRooted = new JLabel("Trees rooted:");
		lblTreeRooted.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblTreeRooted, "flowx,cell 0 5 2 1,alignx right");

		rdbtnTreeRootedYes = new JRadioButton("Yes");
		rdbtnTreeRootedYes.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnTreeRootedYes.setHorizontalAlignment(SwingConstants.LEFT);
		rdbtnTreeRootedYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				TreeRootedToggle(true);
			}
		});
		rdbtnTreeRootedYes.setSelected(false);
		panel.add(rdbtnTreeRootedYes, "cell 2 5");

		rdbtnTreeRootedNo = new JRadioButton("No");
		rdbtnTreeRootedNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnTreeRootedNo.setHorizontalAlignment(SwingConstants.LEFT);
		rdbtnTreeRootedNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				TreeRootedToggle(false);
			}
		});
		rdbtnTreeRootedNo.setSelected(true);
		panel.add(rdbtnTreeRootedNo, "cell 2 5");
		
		lblOutputKind = new JLabel("Output kind:");
		lblOutputKind.setHorizontalAlignment(SwingConstants.RIGHT);
		lblOutputKind.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblOutputKind, "flowx,cell 0 6 2 1,alignx right");
		
		cmbxOutputKind = new JComboBox();
		cmbxOutputKind.setModel(new DefaultComboBoxModel(new String[] {"One pair per line, verbose", "One pair per line, sparse"}));
		cmbxOutputKind.setSelectedIndex(0);
		panel.add(cmbxOutputKind, "cell 2 6, growx");

		lblPrintInd = new JLabel("Display progress:");
		lblPrintInd.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblPrintInd, "flowx,cell 0 7 2 1,alignx right");

		rdbtnPrintIndYes = new JRadioButton("Yes");
		rdbtnPrintIndYes.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnPrintIndYes.setHorizontalAlignment(SwingConstants.LEFT);
		rdbtnPrintIndYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				PrintIndToggle(true);
			}
		});
		rdbtnPrintIndYes.setSelected(true);
		panel.add(rdbtnPrintIndYes, "cell 2 7");

		rdbtnPrintIndNo = new JRadioButton("No");
		rdbtnPrintIndNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnPrintIndNo.setHorizontalAlignment(SwingConstants.LEFT);
		rdbtnPrintIndNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				PrintIndToggle(false);
			}
		});
		rdbtnPrintIndNo.setSelected(false);
		panel.add(rdbtnPrintIndNo, "cell 2 7");

		btnExecute = new JButton("Execute");
		btnExecute.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				inputvals = new TreeDistData();
				inputvals.intree = txtInputTree.getText();
				inputvals.intree2 = txtInputTree2.getText();
				inputvals.outfile = txtOutputFile.getText();
				inputvals.outfileopt = "w";
				
				if (cmbxDistType.getSelectedIndex() == 0)
				{
					inputvals.DistType = "BSD";
				}
				else if (cmbxDistType.getSelectedIndex() == 1)
				{
					inputvals.DistType = "SYMMETRIC";
				}
				else //cmbxDistType.getSelectedIndex() == 2)
				{
					inputvals.DistType = "RF";
				}
				
				if (cmbxDistKind.getSelectedIndex() == 0)
				{
					inputvals.DistKind = "ADJACENT_PAIRS";
					if(cmbxOutputKind.getSelectedIndex() == 0)
					{
						inputvals.OutputKind = "VERBOSE";
					}
					else //if(cmbxOutputKind.getSelectedIndex() == 1)
					{
						inputvals.OutputKind = "SPARSE";	
					}
				}
				else if (cmbxDistKind.getSelectedIndex() == 1)
				{
					inputvals.DistKind = "ALL_IN_FIRST";
					if(cmbxOutputKind.getSelectedIndex() == 0)
					{
						inputvals.OutputKind = "FULL_MATRIX";
					}
					else if(cmbxOutputKind.getSelectedIndex() == 1)
					{
						inputvals.OutputKind = "VERBOSE";
					}
					else // if(cmbxOutputKind.getSelectedIndex() == 2)
					{
						inputvals.OutputKind = "SPARSE";	
					}
				}
				else if (cmbxDistKind.getSelectedIndex() == 2)
				{
					inputvals.DistKind = "CORR_IN_1_AND_2";
					if(cmbxOutputKind.getSelectedIndex() == 0)
					{
						inputvals.OutputKind = "VERBOSE";
					}
					else//if(cmbxOutputKind.getSelectedIndex() == 1)
					{
						inputvals.OutputKind = "SPARSE";	
					}
				}
				else //cmbxDistKind.getSelectedIndex() == 3)
				{
					inputvals.DistKind = "ALL_IN_1_AND_2";
					if(cmbxOutputKind.getSelectedIndex() == 0)
					{
						inputvals.OutputKind = "FULL_MATRIX";
					}
					else if(cmbxOutputKind.getSelectedIndex() == 1)
					{
						inputvals.OutputKind = "VERBOSE";
					}
					else // if(cmbxOutputKind.getSelectedIndex() == 2)
					{
						inputvals.OutputKind = "SPARSE";	
					}
				}
				
				inputvals.Rooted = rdbtnTreeRootedYes.isSelected();
				inputvals.PrintInd = rdbtnPrintIndYes.isSelected();
				
				btnExecute.setEnabled(false);	
				String title = "Treedist Progress";
				String curdir = System.getProperty("user.dir");
				curdir += "/progress.txt";
				File fl = new File(curdir);
				fl.delete();
				inTitle = title;
				inCurdir = curdir;
	
				if (checkInputVals())
				{
					
			  	    Thread treeDistThread = new Thread() {
						public void run() {
							runTreeDistThreads();
						}
			  	    };
			  	    treeDistThread.start();
				}
				btnExecute.setEnabled(true);
			}
		});
		btnExecute.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(btnExecute, "cell 2 8,alignx center");

		btnQuit = new JButton("Quit");
		btnQuit.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				if(phylipCall)
				{
					frmTreeDistControls.dispose();
					return;
				}
				else
				{
					System.exit(0);
				}
			}
		});
		btnQuit.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(btnQuit, "cell 2 8");
	}
	
	public boolean checkInputVals(){
		
		// check files
	    TestFileNames test = new TestFileNames();
		
		if (!test.FileAvailable(inputvals.intree, "Input Tree"))
		{
			return false;
		}
		
		if (inputvals.DistType == "twotree")
		{
			if (!test.FileAvailable(inputvals.intree2, "Input Tree 2"))
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
		
		return true;
	}
	protected void runTreeDistThreads() {
	   	try
	   	{
    		// see if library exists
    		Native.loadLibrary("treedist", TreeDist.class);
		}
    	catch(UnsatisfiedLinkError e)
    	{
     		new TestFileNames().LibraryMissing("TreeDist");
    		return;
    	}
		try 
		{
	  	    Thread treeDistRunThread = new Thread() {
		  	      public void run() {
		  	    	  
		  			// at this point we hook into the C code
		  			TreeDist treedist = (TreeDist) Native.loadLibrary("treedist", TreeDist.class);
		  			treedist.treedist(
		  				inputvals.intree,
		  				inputvals.intree2,
		  				inputvals.outfile,
		  				inputvals.outfileopt,
		  				inputvals.DistType,
		  				inputvals.DistKind,
		  				inputvals.Rooted,
		  				inputvals.OutputKind,
		  				inputvals.PrintInd);
			  	    };
	  	    };
	  	    treeDistRunThread.start();

	  	    if (inputvals.PrintInd)
	  	    {
		  	  	do
		  	  	{
					Thread.sleep(1000);
					updateProgress();
				} while (treeDistRunThread.isAlive());
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
