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

public class DnaInvarUserInterface {
   public interface DnaInvar extends Library {
        public void dnainvar(
        		String infile,
        		String outfile,
        		String outfileopt,
        		String weightfile,
        		boolean SitesWeighted,
        		boolean AnalyzeMult,
        		int NumMult,
        		boolean MultDataSet,
        		boolean InputSeq,
        		boolean PrintData,
        		boolean DotDiff,
        		boolean PrintInd,
        		boolean PrintInvar,
        		boolean PrintPats);
    }

	public class DnaInvarData {
		String infile;
		String outfile;
		String outfileopt;
		String weightfile;
		boolean SitesWeighted;
		boolean AnalyzeMult;
		int NumMult;
		boolean MultDataSet;
		boolean InputSeq;
		boolean PrintData;
		boolean DotDiff;
		boolean PrintInd;
		boolean PrintInvar;
		boolean PrintPats;
	}

	private DnaInvarData inputvals;
	private String inTitle;
	private String inCurdir;
	private DisplayProgress dp;
	private String filedir;
	private boolean ExplicitWgts;
	private boolean phylipCall;

	private JFrame frmDnaInvarControls;
	private JLabel lblSitesWeight;
	private JRadioButton rdbtnSitesYes;
	private JRadioButton rdbtnSitesNo;
	private JLabel lblMultData;
	private JRadioButton rdbtnMultDataYes;
	private JRadioButton rdbtnMultDataNo;
	private JLabel lblMultDataWeight;
	private JRadioButton rdbtnDataSet;
	private JRadioButton rdbtnWeights;
	private JLabel lblInputSeq;
	private JRadioButton rdbtnInputSeqYes;
	private JRadioButton rdbtnInputSeqNo;
	private JLabel lblPrintData;
	private JRadioButton rdbtnPrintDataYes;
	private JRadioButton rdbtnPrintDataNo;
	private JLabel lblPrintInd;
	private JRadioButton rdbtnPrintIndYes;
	private JRadioButton rdbtnPrintIndNo;
	private JLabel lblDotDiff;
	private JRadioButton rdbtnDotDiffYes;
	private JRadioButton rdbtnDotDiffNo;
	private JLabel lblPrintPats;
	private JRadioButton rdbtnPrintPatsYes;
	private JRadioButton rdbtnPrintPatsNo;
	private JLabel lblPrintInvar;
	private JRadioButton rdbtnPrintInvarYes;
	private JRadioButton rdbtnPrintInvarNo;
	private JTextField txtInputFile;
	private JButton btnInputFile;
	private JButton btnExecute;
	private JButton btnQuit;
	private JTextField txtOutputFile;
	private JButton btnOutputFile;
	private JButton btnWeightFile;
	private JTextField txtWeightFile;
	private JLabel lblNSets;
	private JTextField txtNSets;
	
	private JScrollPane scrollPane;
	private JPanel panel;

	/**
	 * Launch the application
	 * 	public static void main(String[] args) {.
	 */
	public static void main(final String[] args) {
		EventQueue.invokeLater(new Runnable() {
			public void run() {
				try {
					DnaInvarUserInterface window = new DnaInvarUserInterface(args);
					window.frmDnaInvarControls.setVisible(true);
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

		int option = fileChooser.showOpenDialog(frmDnaInvarControls.getRootPane());
		if (option == JFileChooser.APPROVE_OPTION) {
			File selectedFile = fileChooser.getSelectedFile();
			filedir = fileChooser.getCurrentDirectory().getAbsolutePath();
			file.setText(selectedFile.getPath());
		}
	}

	protected void SiteToggle(boolean isSite) {
		if (isSite) {
			rdbtnSitesYes.setSelected(true);
			rdbtnSitesNo.setSelected(false);
			btnWeightFile.setEnabled(true);
			txtWeightFile.setEnabled(true);
			ExplicitWgts = true;
		} else {
			rdbtnSitesYes.setSelected(false);
			rdbtnSitesNo.setSelected(true);
			btnWeightFile.setEnabled(false);
			txtWeightFile.setEnabled(false);
			ExplicitWgts = false;
		}
	}

	protected void MultToggle(boolean isMult) {
		if (isMult) {
			rdbtnMultDataYes.setSelected(true);
			rdbtnMultDataNo.setSelected(false);
			lblMultDataWeight.setEnabled(true);
			rdbtnDataSet.setEnabled(true);
			rdbtnWeights.setEnabled(true);
			lblNSets.setEnabled(true);
			txtNSets.setEnabled(true);
			lblInputSeq.setEnabled(true);
			rdbtnInputSeqYes.setEnabled(true);
			rdbtnInputSeqNo.setEnabled(true);
			if (rdbtnWeights.isSelected())
			{
				if(!ExplicitWgts)
				{
					btnWeightFile.setEnabled(true);
					txtWeightFile.setEnabled(true);	
					rdbtnSitesYes.setSelected(true);
					rdbtnSitesNo.setSelected(false);
				}
			}
		} else {
			rdbtnMultDataYes.setSelected(false);
			rdbtnMultDataNo.setSelected(true);
			lblMultDataWeight.setEnabled(false);
			rdbtnDataSet.setEnabled(false);
			rdbtnWeights.setEnabled(false);
			lblNSets.setEnabled(false);
			txtNSets.setEnabled(false);
			lblInputSeq.setEnabled(false);
			rdbtnInputSeqYes.setEnabled(false);
			rdbtnInputSeqNo.setEnabled(false);
			if(ExplicitWgts)
			{
				btnWeightFile.setEnabled(true);
				txtWeightFile.setEnabled(true);	
				rdbtnSitesYes.setSelected(true);
				rdbtnSitesNo.setSelected(false);
			}
			else
			{
				btnWeightFile.setEnabled(false);
				txtWeightFile.setEnabled(false);	
				rdbtnSitesYes.setSelected(false);
				rdbtnSitesNo.setSelected(true);				
			}
		}
	}

	protected void DataWeightToggle(boolean isyes) {
		if (isyes) {
			rdbtnDataSet.setSelected(true);
			rdbtnWeights.setSelected(false);
			if(!ExplicitWgts)
			{
				btnWeightFile.setEnabled(false);
				txtWeightFile.setEnabled(false);
				rdbtnSitesYes.setSelected(false);
				rdbtnSitesNo.setSelected(true);
			}
		} else {
			rdbtnDataSet.setSelected(false);
			rdbtnWeights.setSelected(true);
			btnWeightFile.setEnabled(true);
			txtWeightFile.setEnabled(true);
			rdbtnSitesYes.setSelected(true);
			rdbtnSitesNo.setSelected(false);
		}
	}

	protected void InputSeqToggle(boolean isInput) {
		if (isInput) {
			rdbtnInputSeqYes.setSelected(true);
			rdbtnInputSeqNo.setSelected(false);
		} else {
			rdbtnInputSeqYes.setSelected(false);
			rdbtnInputSeqNo.setSelected(true);
		}
	}

	protected void PrintDataToggle(boolean isPrintData) {
		if (isPrintData) {
			rdbtnPrintDataYes.setSelected(true);
			rdbtnPrintDataNo.setSelected(false);
			lblDotDiff.setEnabled(true);
			rdbtnDotDiffYes.setEnabled(true);
			rdbtnDotDiffNo.setEnabled(true);
		} else {
			rdbtnPrintDataYes.setSelected(false);
			rdbtnPrintDataNo.setSelected(true);
			lblDotDiff.setEnabled(false);
			rdbtnDotDiffYes.setEnabled(false);
			rdbtnDotDiffNo.setEnabled(false);
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

	protected void DotDiffToggle(boolean isDot) {
		if (isDot) {
			rdbtnDotDiffYes.setSelected(true);
			rdbtnDotDiffNo.setSelected(false);
		} else {
			rdbtnDotDiffYes.setSelected(false);
			rdbtnDotDiffNo.setSelected(true);
		}
	}

	protected void PrintSeqToggle(boolean isSeq) {
		if (isSeq) {
			rdbtnPrintPatsYes.setSelected(true);
			rdbtnPrintPatsNo.setSelected(false);
		} else {
			rdbtnPrintPatsYes.setSelected(false);
			rdbtnPrintPatsNo.setSelected(true);
		}
	}
	
	/**
	 * Create the application.
	 */
	public DnaInvarUserInterface(String[] args) {
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
		ExplicitWgts = false;

		frmDnaInvarControls = new JFrame();
		frmDnaInvarControls.setBackground(new Color(204, 255, 255));
		frmDnaInvarControls.setTitle("Dnainvar");
		frmDnaInvarControls.setBounds(100, 100, 640, 470);
		frmDnaInvarControls.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frmDnaInvarControls.setPreferredSize(new Dimension(frmDnaInvarControls.getBounds().width, frmDnaInvarControls.getBounds().height));
		
		scrollPane = new JScrollPane();
		scrollPane.setPreferredSize(frmDnaInvarControls.getPreferredSize());
		frmDnaInvarControls.getContentPane().add(scrollPane);
		
		panel = new JPanel();
		panel.setPreferredSize(frmDnaInvarControls.getPreferredSize());
		scrollPane.setViewportView(panel);
		panel.setLayout(new MigLayout("", "[pref!,grow][10.00,grow][pref!,grow]", "[][][][]"));

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
		panel.add(txtInputFile, "cell 1 0 2 1,growx");
		
		btnWeightFile = new JButton("Weights File");
		btnWeightFile.setEnabled(false);
		btnWeightFile.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				ChooseFile(txtWeightFile);
			}
		});
		panel.add(btnWeightFile, "cell 0 1,growx");
		
		txtWeightFile = new JTextField();
		txtWeightFile.setEnabled(false);
		txtWeightFile.setFont(new Font("Arial", Font.PLAIN, 13));
		txtWeightFile.setText("weightfile");
		panel.add(txtWeightFile, "cell 1 1 2 1,growx");
		
		btnOutputFile = new JButton("Output File");
		btnOutputFile.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				ChooseFile(txtOutputFile);
			}
		});
		panel.add(btnOutputFile, "cell 0 2,growx");
		
		txtOutputFile = new JTextField();
		txtOutputFile.setFont(new Font("Arial", Font.PLAIN, 13));
		txtOutputFile.setText("outfile");
		panel.add(txtOutputFile, "cell 1 2 2 1,growx");

		lblSitesWeight = new JLabel("Sites weighted:");
		lblSitesWeight.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblSitesWeight, "flowx,cell 0 3 2 1,alignx right");

		rdbtnSitesYes = new JRadioButton("Yes");
		rdbtnSitesYes.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnSitesYes.setHorizontalAlignment(SwingConstants.LEFT);
		rdbtnSitesYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				SiteToggle(true);
			}
		});
		panel.add(rdbtnSitesYes, "cell 2 3");

		rdbtnSitesNo = new JRadioButton("No");
		rdbtnSitesNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnSitesNo.setHorizontalAlignment(SwingConstants.LEFT);
		rdbtnSitesNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				SiteToggle(false);
			}
		});
		rdbtnSitesNo.setSelected(true);
		panel.add(rdbtnSitesNo, "cell 2 3");

		lblMultData = new JLabel("Analyze multiple data sets:");
		lblMultData.setFont(new Font("Arial", Font.BOLD, 13));
		lblMultData.setHorizontalAlignment(SwingConstants.RIGHT);
		panel.add(lblMultData, "flowx,cell 0 4 2 1,alignx right");

		rdbtnMultDataYes = new JRadioButton("Yes");
		rdbtnMultDataYes.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnMultDataYes.setHorizontalAlignment(SwingConstants.LEFT);
		rdbtnMultDataYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				MultToggle(true);
			}
		});
		rdbtnMultDataYes.setSelected(false);
		panel.add(rdbtnMultDataYes, "cell 2 4");

		rdbtnMultDataNo = new JRadioButton("No");
		rdbtnMultDataNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnMultDataNo.setHorizontalAlignment(SwingConstants.LEFT);
		rdbtnMultDataNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				MultToggle(false);
			}
		});
		rdbtnMultDataNo.setSelected(true);
		panel.add(rdbtnMultDataNo, "cell 2 4");

		lblMultDataWeight = new JLabel(
				"Multiple data sets or multiple weights:");
		lblMultDataWeight.setFont(new Font("Arial", Font.BOLD, 13));
		lblMultDataWeight.setHorizontalAlignment(SwingConstants.RIGHT);
		lblMultDataWeight.setEnabled(false);
		panel.add(lblMultDataWeight, "flowx,cell 0 5 2 1,alignx right");

		rdbtnDataSet = new JRadioButton("Data sets");
		rdbtnDataSet.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnDataSet.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				DataWeightToggle(true);
			}
		});
		rdbtnDataSet.setSelected(true);
		rdbtnDataSet.setEnabled(false);
		panel.add(rdbtnDataSet, "cell 2 5");

		rdbtnWeights = new JRadioButton("Weights");
		rdbtnWeights.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnWeights.setHorizontalAlignment(SwingConstants.LEFT);
		rdbtnWeights.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				DataWeightToggle(false);
			}
		});
		rdbtnWeights.setSelected(false);
		rdbtnWeights.setEnabled(false);
		panel.add(rdbtnWeights, "cell 2 5");
		
		lblNSets = new JLabel("Number:\r\n");
		lblNSets.setHorizontalAlignment(SwingConstants.RIGHT);
		lblNSets.setFont(new Font("Arial", Font.BOLD, 13));
		lblNSets.setEnabled(false);
		panel.add(lblNSets, "cell 2 6");
		
		txtNSets = new JTextField();
		txtNSets.setText("1");
		txtNSets.setFont(new Font("Arial", Font.PLAIN, 13));
		txtNSets.setEnabled(false);
		txtNSets.setColumns(6);
		panel.add(txtNSets, "cell 2 6");

		lblInputSeq = new JLabel("Input sequences:");
		lblInputSeq.setEnabled(false);
		lblInputSeq.setFont(new Font("Arial", Font.BOLD, 13));
		lblInputSeq.setHorizontalAlignment(SwingConstants.RIGHT);
		panel.add(lblInputSeq, "flowx,cell 0 7 2 1,alignx right");

		rdbtnInputSeqYes = new JRadioButton("Interleaved");
		rdbtnInputSeqYes.setEnabled(false);
		rdbtnInputSeqYes.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnInputSeqYes.setHorizontalAlignment(SwingConstants.LEFT);
		rdbtnInputSeqYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				InputSeqToggle(true);
			}
		});
		rdbtnInputSeqYes.setSelected(true);
		panel.add(rdbtnInputSeqYes, "cell 2 7");

		rdbtnInputSeqNo = new JRadioButton("Sequential");
		rdbtnInputSeqNo.setEnabled(false);
		rdbtnInputSeqNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnInputSeqNo.setHorizontalAlignment(SwingConstants.LEFT);
		rdbtnInputSeqNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				InputSeqToggle(false);
			}
		});
		rdbtnInputSeqNo.setSelected(false);
		panel.add(rdbtnInputSeqNo, "cell 2 7");

		lblPrintData = new JLabel("Print out the data at start of run:");
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
		rdbtnPrintDataYes.setSelected(false);
		panel.add(rdbtnPrintDataYes, "cell 2 8");

		rdbtnPrintDataNo = new JRadioButton("No");
		rdbtnPrintDataNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnPrintDataNo.setHorizontalAlignment(SwingConstants.LEFT);
		rdbtnPrintDataNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				PrintDataToggle(false);
			}
		});
		rdbtnPrintDataNo.setSelected(true);
		panel.add(rdbtnPrintDataNo, "cell 2 8");

		lblDotDiff = new JLabel("Use dot-differencing to display:");
		lblDotDiff.setEnabled(false);
		lblDotDiff.setFont(new Font("Arial", Font.BOLD, 13));
		lblDotDiff.setHorizontalAlignment(SwingConstants.RIGHT);
		panel.add(lblDotDiff, "flowx,cell 0 9 2 1,alignx right");

		rdbtnDotDiffYes = new JRadioButton("Yes");
		rdbtnDotDiffYes.setEnabled(false);
		rdbtnDotDiffYes.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnDotDiffYes.setHorizontalAlignment(SwingConstants.LEFT);
		rdbtnDotDiffYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				DotDiffToggle(true);
			}
		});
		rdbtnDotDiffYes.setSelected(true);
		panel.add(rdbtnDotDiffYes, "cell 2 9");

		rdbtnDotDiffNo = new JRadioButton("No");
		rdbtnDotDiffNo.setEnabled(false);
		rdbtnDotDiffNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnDotDiffNo.setHorizontalAlignment(SwingConstants.LEFT);
		rdbtnDotDiffNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				DotDiffToggle(false);
			}
		});
		panel.add(rdbtnDotDiffNo, "cell 2 9");

		lblPrintPats = new JLabel("Print counts of patterns:");
		lblPrintPats.setFont(new Font("Arial", Font.BOLD, 13));
		lblPrintPats.setHorizontalAlignment(SwingConstants.RIGHT);
		panel.add(lblPrintPats, "flowx,cell 0 10 2 1,alignx right");

		rdbtnPrintPatsYes = new JRadioButton("Yes");
		rdbtnPrintPatsYes.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnPrintPatsYes.setHorizontalAlignment(SwingConstants.LEFT);
		rdbtnPrintPatsYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				PrintSeqToggle(true);
			}
		});
		rdbtnPrintPatsYes.setSelected(true);
		panel.add(rdbtnPrintPatsYes, "cell 2 10");

		rdbtnPrintPatsNo = new JRadioButton("No");
		rdbtnPrintPatsNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnPrintPatsNo.setHorizontalAlignment(SwingConstants.LEFT);
		rdbtnPrintPatsNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				PrintSeqToggle(false);
			}
		});
		panel.add(rdbtnPrintPatsNo, "cell 2 10");

		lblPrintInvar = new JLabel("Print the invariants:");
		lblPrintInvar.setFont(new Font("Arial", Font.BOLD, 13));
		lblPrintInvar.setHorizontalAlignment(SwingConstants.RIGHT);
		panel.add(lblPrintInvar, "flowx,cell 0 11 2 1,alignx right");

		rdbtnPrintInvarYes = new JRadioButton("Yes");
		rdbtnPrintInvarYes.setSelected(true);
		rdbtnPrintInvarYes.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnPrintInvarYes.setHorizontalAlignment(SwingConstants.LEFT);
		panel.add(rdbtnPrintInvarYes, "cell 2 11");

		rdbtnPrintInvarNo = new JRadioButton("No");
		rdbtnPrintInvarNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnPrintInvarNo.setHorizontalAlignment(SwingConstants.LEFT);
		panel.add(rdbtnPrintInvarNo, "cell 2 11");

		lblPrintInd = new JLabel("Display progress:");
		lblPrintInd.setFont(new Font("Arial", Font.BOLD, 13));
		lblPrintInd.setHorizontalAlignment(SwingConstants.RIGHT);
		panel.add(lblPrintInd, "flowx,cell 0 12 2 1,alignx right");

		rdbtnPrintIndYes = new JRadioButton("Yes");
		rdbtnPrintIndYes.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnPrintIndYes.setHorizontalAlignment(SwingConstants.LEFT);
		rdbtnPrintIndYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				PrintIndToggle(true);
			}
		});
		rdbtnPrintIndYes.setSelected(true);
		panel.add(rdbtnPrintIndYes, "cell 2 12");

		rdbtnPrintIndNo = new JRadioButton("No");
		rdbtnPrintIndNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnPrintIndNo.setHorizontalAlignment(SwingConstants.LEFT);
		rdbtnPrintIndNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				PrintIndToggle(false);
			}
		});
		rdbtnPrintIndNo.setSelected(false);
		panel.add(rdbtnPrintIndNo, "cell 2 12");

		btnExecute = new JButton("Execute");
		btnExecute.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				inputvals = new DnaInvarData();
				inputvals.infile = txtInputFile.getText();
				inputvals.outfile = txtOutputFile.getText();
				inputvals.outfileopt = "w";
				inputvals.weightfile = txtWeightFile.getText();
				inputvals.SitesWeighted = rdbtnSitesYes.isSelected();
				inputvals.AnalyzeMult = rdbtnMultDataYes.isSelected();
				inputvals.MultDataSet = rdbtnDataSet.isSelected();
				inputvals.NumMult = Integer.parseInt(txtNSets.getText());
				inputvals.InputSeq = rdbtnInputSeqYes.isSelected();
				inputvals.PrintData = rdbtnPrintDataYes.isSelected();
				inputvals.DotDiff = rdbtnDotDiffYes.isSelected();
				inputvals.PrintInd = rdbtnPrintIndYes.isSelected();
				inputvals.PrintPats = rdbtnPrintPatsYes.isSelected();
				inputvals.PrintInvar = rdbtnPrintInvarYes.isSelected();
				
				btnExecute.setEnabled(false);
				String title = "Dnainvar Progress";
				String curdir = System.getProperty("user.dir");
				curdir += "/progress.txt";
				File fl = new File(curdir);
				fl.delete();
				inTitle = title;
				inCurdir = curdir;
	
				if (checkInputVals())
				{
					
			  	    Thread dnaInvarThread = new Thread() {
						public void run() {
							runDnaInvarThreads();
						}
			  	    };
			  	  dnaInvarThread.start();
				}
				btnExecute.setEnabled(true);
			}
		});
		btnExecute.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(btnExecute, "cell 2 13");

		btnQuit = new JButton("Quit");
		btnQuit.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				if(phylipCall)
				{
					frmDnaInvarControls.dispose();
					return;
				}
				else
				{
					System.exit(0);
				}
			}
		});
		btnQuit.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(btnQuit, "cell 2 13");

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
		
		String opt = test.FileAlreadyExists(inputvals.outfile, "Outfile");
		if (opt == "q")
		{
			return false;
		}
		else
		{
			inputvals.outfileopt = opt;
		}
		
		if (inputvals.SitesWeighted)
		{
			if (!test.FileAvailable(inputvals.weightfile, "Weights"))
			{
				return false;
			}
		}
		return true;
	}
	
	protected void runDnaInvarThreads() {
    	try
    	{
    		// see if library exists
    		Native.loadLibrary("dnainvar", DnaInvar.class);
		}
    	catch(UnsatisfiedLinkError e)
    	{
     		new TestFileNames().LibraryMissing("DnaInvar");
    		return;
    	}
		try 
		{
	  	    Thread dnaInvarRunThread = new Thread() {
		  	      public void run() {
		  	    	  
		  			// at this point we hook into the C code
		  			DnaInvar dnainvar = (DnaInvar) Native.loadLibrary("dnainvar", DnaInvar.class);
		  			dnainvar.dnainvar(
		  	     		inputvals.infile,
		  	     		inputvals.outfile,
		  	     		inputvals.outfileopt,
		  	     		inputvals.weightfile,
		  	     		inputvals.SitesWeighted,
		  	     		inputvals.AnalyzeMult,
		  	    		inputvals.NumMult,
		  	    		inputvals.MultDataSet,
		  	     		inputvals.InputSeq,
		  	     		inputvals.PrintData,
		  	     		inputvals.DotDiff,
		  	     		inputvals.PrintInd,
		  	     		inputvals.PrintInvar,
		  	     		inputvals.PrintPats);
		  	    };
	  	    };
	  	    dnaInvarRunThread.start();

	  	    if (inputvals.PrintInd)
	  	    {
		  	  	do
		  	  	{
					Thread.sleep(1000);
					updateProgress();
				} while (dnaInvarRunThread.isAlive());
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
