package phylip;
import java.awt.EventQueue;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JRadioButton;
import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;
import java.io.File;
import javax.swing.JFileChooser;
import javax.swing.JOptionPane;
import javax.swing.JTextField;
import javax.swing.JButton;
import javax.swing.SwingConstants;
import javax.swing.SwingUtilities;

import com.sun.jna.Library;
import com.sun.jna.Native;

import utilities.DisplayProgress;
import utilities.TestFileNames;

import java.awt.Font;
import java.awt.Color;

import net.miginfocom.swing.MigLayout;
import javax.swing.JScrollPane;
import java.awt.Dimension;
import javax.swing.JPanel;

public class DollopUserInterface {
   public interface Dollop extends Library {
        public void dollop(
           		String infile,
           		String intree,
           		String weightfile,
           		String ancfile,
           		String outfile,
           		String outfileopt,
           		String outtree,
           		String outtreeopt,
           		boolean BestTree,
           		boolean UseDollo,
           		int TreeMax,
           		boolean InputOrder,
           		int RandNum,
           		int NumJumble,
           		boolean ThreshDollop,
           		double ThreshVal,
           		boolean UseAncStates,
           		boolean SitesWeighted,
           		boolean AnalyzeMult,
           		boolean MultDataset,
           		int NumMult,
           		boolean PrintData,
           		boolean PrintInd,
           		boolean PrintTree,
           		boolean PrintSteps,
           		boolean PrintSeq,
           		boolean WriteTree);
        }
	
	public class DollopData {
		String infile;
		String intree;
		String weightfile;
		String ancfile;
		String outfile;
		String outfileopt;
		String outtree;
		String outtreeopt;
		boolean BestTree;
		boolean UseDollo;
		int TreeMax;
		boolean InputOrder;
		int RandNum;
		int NumJumble;
		boolean ThreshDollop;
		double ThreshVal;
		boolean UseAncStates;
		boolean SitesWeighted;
		boolean AnalyzeMult;
		boolean MultDataset;
		int NumMult;
		boolean PrintData;
		boolean PrintInd;
		boolean PrintTree;
		boolean PrintSteps;
		boolean PrintSeq;
		boolean WriteTree;
	}

	private DollopData inputvals;
	private String inTitle;
	private String inCurdir;
	private DisplayProgress dp;
	private String filedir;
	private boolean ExplicitWgts;
	private boolean phylipCall;

	private JFrame frmDollopControls;
	private JTextField txtInputFile;
	private JButton btnInputFile;
	private JTextField txtOutputFile;
	private JButton btnOutputFile;
	private JButton btnWeightFile;
	private JTextField txtWeightFile;
	private JButton btnAncFile;
	private JTextField txtAncFile;
	private JButton btnOutputTree;
	private JTextField txtOutputTree;
	private JLabel lblSearchTree;
	private JRadioButton rbtnSearchYes;
	private JRadioButton rbtnSearchNo;
	private JLabel lblMethod;
	private JRadioButton rbtnDollo;
	private JRadioButton rbtnPoly;
	private JLabel lblRandOrder;
	private JRadioButton rbtnRandOrderYes;
	private JRadioButton rbtnRandOrderNo;
	private JLabel lblRandomSeed;
	private JTextField txtRandomSeed;
	private JLabel lblNumberJumble;
	private JTextField txtNumberJumble;
	private JLabel lblUseThreshold;
	private JRadioButton rbtnThresholdYes;
	private JRadioButton rbtnThresholdNo;
	private JLabel lblThresholdValue;
	private JTextField txtThresholdValue;
	private JLabel lblSitesWeight;
	private JRadioButton rbtnSitesYes;
	private JRadioButton rbtnSitesNo;
	private JLabel lblMultData;
	private JRadioButton rbtnMultDataYes;
	private JRadioButton rbtnMultDataNo;
	private JLabel lblMultDataWeight;
	private JRadioButton rbtnDataSet;
	private JRadioButton rbtnWeights;
	private JLabel lblPrintData;
	private JRadioButton rbtnPrintDataYes;
	private JRadioButton rbtnPrintDataNo;
	private JLabel lblPrintInd;
	private JRadioButton rbtnPrintIndYes;
	private JRadioButton rbtnPrintIndNo;
	private JLabel lblPrintTree;
	private JRadioButton rbtnPrintTreeYes;
	private JRadioButton rbtnPrintTreeNo;
	private JLabel lblPrintSteps;
	private JRadioButton rbtnPrintStepYes;
	private JRadioButton rbtnPrintStepNo;
	private JLabel lblPrintChar;
	private JRadioButton rbtnPrintCharYes;
	private JRadioButton rbtnPrintCharNo;
	private JLabel lblWriteTree;
	private JRadioButton rbtnWriteTreeYes;
	private JRadioButton rbtnWriteTreeNo;
	private JButton btnExecute;
	private JButton btnQuit;
	private JLabel lblNSets;
	private JTextField txtNSets;
	private JLabel lblRandOdd;
	private JLabel lblAncStates;
	private JRadioButton rbtnAncStatesYes;
	private JRadioButton rbtnAncStatesNo;
	private JButton btnInTree;
	private JTextField txtInTree;
	private JLabel lblNumberOfTrees;
	private JTextField txtNumberOfTrees;
	
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
					DollopUserInterface window = new DollopUserInterface(args);
					window.frmDollopControls.setVisible(true);
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

		int option = fileChooser.showOpenDialog(frmDollopControls.getRootPane());
		if (option == JFileChooser.APPROVE_OPTION) {
			File selectedFile = fileChooser.getSelectedFile();
			filedir = fileChooser.getCurrentDirectory().getAbsolutePath();
			file.setText(selectedFile.getPath());
		}
	}

	protected void SearchToggle(boolean isSearch) {
		if (isSearch) {
			rbtnSearchYes.setSelected(true);
			rbtnSearchNo.setSelected(false);
			btnInTree.setEnabled(false);
			txtInTree.setEnabled(false);
			lblRandOrder.setEnabled(true);
			rbtnRandOrderYes.setEnabled(true);
			rbtnRandOrderNo.setEnabled(true);
			RandToggle(rbtnRandOrderNo.isSelected());
		} else {
			rbtnSearchYes.setSelected(false);
			rbtnSearchNo.setSelected(true);
			btnInTree.setEnabled(true);
			txtInTree.setEnabled(true);
			lblRandOrder.setEnabled(false);
			rbtnRandOrderYes.setEnabled(false);
			rbtnRandOrderNo.setEnabled(false);
			lblRandomSeed.setEnabled(false);
			txtRandomSeed.setEnabled(false);
			lblNumberJumble.setEnabled(false);
			txtNumberJumble.setEnabled(false);
			lblRandOdd.setEnabled(false);
		}
	}
	
	protected void RandToggle(boolean isRand) {
		if (isRand) {
			rbtnRandOrderYes.setSelected(false);
			rbtnRandOrderNo.setSelected(true);
			lblRandomSeed.setEnabled(false);
			txtRandomSeed.setEnabled(false);
			lblNumberJumble.setEnabled(false);
			txtNumberJumble.setEnabled(false);
			lblRandOdd.setEnabled(false);
		} else {
			rbtnRandOrderYes.setSelected(true);
			rbtnRandOrderNo.setSelected(false);
			lblRandomSeed.setEnabled(true);
			txtRandomSeed.setEnabled(true);
			lblNumberJumble.setEnabled(true);
			txtNumberJumble.setEnabled(true);
			lblRandOdd.setEnabled(true);
		}
	}

	protected void SiteToggle(boolean isSite) {
		if (isSite) {
			rbtnSitesYes.setSelected(true);
			rbtnSitesNo.setSelected(false);
			btnWeightFile.setEnabled(true);
			txtWeightFile.setEnabled(true);
			ExplicitWgts = true;
		} else {
			rbtnSitesYes.setSelected(false);
			rbtnSitesNo.setSelected(true);
			btnWeightFile.setEnabled(false);
			txtWeightFile.setEnabled(false);
			ExplicitWgts = false;
		}
	}

	protected void MultToggle(boolean isMult) {
		if (isMult) {
			rbtnMultDataYes.setSelected(true);
			rbtnMultDataNo.setSelected(false);
			lblMultDataWeight.setEnabled(true);
			rbtnDataSet.setEnabled(true);
			rbtnWeights.setEnabled(true);
			lblNSets.setEnabled(true);
			txtNSets.setEnabled(true);
			if (rbtnWeights.isSelected())
			{
				if(!ExplicitWgts)
				{
					btnWeightFile.setEnabled(true);
					txtWeightFile.setEnabled(true);	
					rbtnSitesYes.setSelected(true);
					rbtnSitesNo.setSelected(false);
				}
			}
		} else {
			rbtnMultDataYes.setSelected(false);
			rbtnMultDataNo.setSelected(true);
			lblMultDataWeight.setEnabled(false);
			rbtnDataSet.setEnabled(false);
			rbtnWeights.setEnabled(false);
			lblNSets.setEnabled(false);
			txtNSets.setEnabled(false);
			if(ExplicitWgts)
			{
				btnWeightFile.setEnabled(true);
				txtWeightFile.setEnabled(true);	
				rbtnSitesYes.setSelected(true);
				rbtnSitesNo.setSelected(false);
			}
			else
			{
				btnWeightFile.setEnabled(false);
				txtWeightFile.setEnabled(false);	
				rbtnSitesYes.setSelected(false);
				rbtnSitesNo.setSelected(true);				
			}
		}
	}

	protected void DataWeightToggle(boolean isyes) {
		if (isyes) {
			rbtnDataSet.setSelected(true);
			rbtnWeights.setSelected(false);
			if(!ExplicitWgts)
			{
				btnWeightFile.setEnabled(false);
				txtWeightFile.setEnabled(false);
				rbtnSitesYes.setSelected(false);
				rbtnSitesNo.setSelected(true);
			}
		} else {
			rbtnDataSet.setSelected(false);
			rbtnWeights.setSelected(true);
			btnWeightFile.setEnabled(true);
			txtWeightFile.setEnabled(true);
			rbtnSitesYes.setSelected(true);
			rbtnSitesNo.setSelected(false);
		}
	}

	protected void PrintDataToggle(boolean isPrintData) {
		if (isPrintData) {
			rbtnPrintDataYes.setSelected(true);
			rbtnPrintDataNo.setSelected(false);
		} else {
			rbtnPrintDataYes.setSelected(false);
			rbtnPrintDataNo.setSelected(true);
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

	protected void PrintTreeToggle(boolean isTree) {
		if (isTree) {
			rbtnPrintTreeYes.setSelected(true);
			rbtnPrintTreeNo.setSelected(false);
		} else {
			rbtnPrintTreeYes.setSelected(false);
			rbtnPrintTreeNo.setSelected(true);
		}
	}

	protected void PrintStepToggle(boolean isStep) {
		if (isStep) {
			rbtnPrintStepYes.setSelected(true);
			rbtnPrintStepNo.setSelected(false);
		} else {
			rbtnPrintStepYes.setSelected(false);
			rbtnPrintStepNo.setSelected(true);
		}
	}

	protected void PrintSeqToggle(boolean isSeq) {
		if (isSeq) {
			rbtnPrintCharYes.setSelected(true);
			rbtnPrintCharNo.setSelected(false);
		} else {
			rbtnPrintCharYes.setSelected(false);
			rbtnPrintCharNo.setSelected(true);
		}
	}

	protected void WriteTreeToggle(boolean isWrite) {
		if (isWrite) {
			rbtnWriteTreeYes.setSelected(true);
			rbtnWriteTreeNo.setSelected(false);
			btnOutputTree.setEnabled(true);
			txtOutputTree.setEnabled(true);
		} else {
			rbtnWriteTreeYes.setSelected(false);
			rbtnWriteTreeNo.setSelected(true);
			btnOutputTree.setEnabled(false);
			txtOutputTree.setEnabled(false);
		}
	}

	protected void ThresholdToggle(boolean isThresh) {
		if (isThresh) {
			rbtnThresholdYes.setSelected(true);
			rbtnThresholdNo.setSelected(false);
			lblThresholdValue.setEnabled(true);
			txtThresholdValue.setEnabled(true);
		} else {
			rbtnThresholdYes.setSelected(false);
			rbtnThresholdNo.setSelected(true);
			lblThresholdValue.setEnabled(false);
			txtThresholdValue.setEnabled(false);
		}
	}
	
	protected void ReadAncestorsToggle(boolean isYes){
		if (isYes)
		{
			rbtnAncStatesYes.setSelected(true);
			rbtnAncStatesNo.setSelected(false);
			btnAncFile.setEnabled(true);
			txtAncFile.setEnabled(true);
		} else {
			rbtnAncStatesYes.setSelected(false);
			rbtnAncStatesNo.setSelected(true);
			btnAncFile.setEnabled(false);
			txtAncFile.setEnabled(false);
		}
	}
	
	protected void MethodToggle(boolean isDollo) {
		if (isDollo) {
			rbtnDollo.setSelected(true);
			rbtnPoly.setSelected(false);
		} else {
			rbtnDollo.setSelected(false);
			rbtnPoly.setSelected(true);
		}
	}



	/**
	 * Create the application.
	 */
	public DollopUserInterface(String[] args) {
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

		frmDollopControls = new JFrame();
		frmDollopControls.setBackground(new Color(204, 255, 255));
		frmDollopControls.setTitle("Dollop");
		frmDollopControls.setBounds(100, 100, 700, 800);
		frmDollopControls.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frmDollopControls.setPreferredSize(new Dimension(frmDollopControls.getBounds().width, frmDollopControls.getBounds().height));
		
		scrollPane = new JScrollPane();
		scrollPane.setPreferredSize(frmDollopControls.getPreferredSize());
		frmDollopControls.getContentPane().add(scrollPane);
		
		panel = new JPanel();
		panel.setPreferredSize(frmDollopControls.getPreferredSize());
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
		panel.add(txtInputFile, "cell 1 0 2 1,growx");
		
		btnInTree = new JButton("Input Tree");
		btnInTree.setFont(new Font("Arial", Font.BOLD, 13));
		btnInTree.setEnabled(false);
		panel.add(btnInTree, "cell 0 1,growx");
		
		txtInTree = new JTextField();
		txtInTree.setText("intree");
		txtInTree.setFont(new Font("Arial", Font.PLAIN, 13));
		txtInTree.setEnabled(false);
		panel.add(txtInTree, "cell 1 1 2 1,growx");
	
		btnWeightFile = new JButton("Weights File");
		btnWeightFile.setEnabled(false);
		btnWeightFile.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				ChooseFile(txtWeightFile);
			}
		});
		btnWeightFile.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(btnWeightFile, "cell 0 2,growx");
		
		txtWeightFile = new JTextField();
		txtWeightFile.setEnabled(false);
		txtWeightFile.setFont(new Font("Arial", Font.PLAIN, 13));
		txtWeightFile.setText("weightfile");
		panel.add(txtWeightFile, "cell 1 2 2 1,growx");
		
		btnAncFile = new JButton("Ancestor File");
		btnAncFile.setFont(new Font("Arial", Font.BOLD, 13));
		btnAncFile.setEnabled(false);
		btnAncFile.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				ChooseFile(txtAncFile);
			}
		});
		panel.add(btnAncFile, "cell 0 3,growx");
		
		txtAncFile = new JTextField();
		txtAncFile.setText("ancestorfile");
		txtAncFile.setFont(new Font("Arial", Font.PLAIN, 13));
		txtAncFile.setEnabled(false);
		panel.add(txtAncFile, "cell 1 3 2 1,growx");
		
		btnOutputFile = new JButton("Output File");
		btnOutputFile.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				ChooseFile(txtOutputFile);
			}
		});
		btnOutputFile.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(btnOutputFile, "cell 0 4,growx");
		
		txtOutputFile = new JTextField();
		txtOutputFile.setFont(new Font("Arial", Font.PLAIN, 13));
		txtOutputFile.setText("outfile");
		panel.add(txtOutputFile, "cell 1 4 2 1,growx");
		
		btnOutputTree = new JButton("Output Tree");
		btnOutputTree.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				ChooseFile(txtOutputTree);
			}
		});
		btnOutputTree.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(btnOutputTree, "cell 0 5,growx");
		
		txtOutputTree = new JTextField();
		txtOutputTree.setText("outtree");
		txtOutputTree.setFont(new Font("Arial", Font.PLAIN, 13));
		panel.add(txtOutputTree, "cell 1 5 2 1,growx");

		lblSearchTree = new JLabel("Search for best tree:");
		lblSearchTree.setFont(new Font("Arial", Font.BOLD, 13));
		lblSearchTree.setHorizontalAlignment(SwingConstants.RIGHT);
		panel.add(lblSearchTree, "flowx,cell 0 6 2 1,alignx right");

		rbtnSearchYes = new JRadioButton("Yes");
		rbtnSearchYes.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnSearchYes.setHorizontalAlignment(SwingConstants.LEFT);
		rbtnSearchYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				SearchToggle(true);
			}
		});
		rbtnSearchYes.setSelected(true);
		panel.add(rbtnSearchYes, "cell 2 6");

		rbtnSearchNo = new JRadioButton("No, use user trees in input file");
		rbtnSearchNo.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnSearchNo.setHorizontalAlignment(SwingConstants.LEFT);
		rbtnSearchNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				SearchToggle(false);
			}
		});
		rbtnSearchNo.setSelected(false);
		panel.add(rbtnSearchNo, "cell 2 6");
		
		lblNumberOfTrees = new JLabel("Number of trees:");
		lblNumberOfTrees.setHorizontalAlignment(SwingConstants.RIGHT);
		lblNumberOfTrees.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblNumberOfTrees, "flowx,cell 0 7 2 1,alignx right");
		
		txtNumberOfTrees = new JTextField();
		txtNumberOfTrees.setText("10000");
		txtNumberOfTrees.setFont(new Font("Arial", Font.PLAIN, 13));
		txtNumberOfTrees.setColumns(6);
		panel.add(txtNumberOfTrees, "cell 2 7");

		lblMethod = new JLabel("Parsimony method:");
		lblMethod.setFont(new Font("Arial", Font.BOLD, 13));
		lblMethod.setHorizontalAlignment(SwingConstants.RIGHT);
		panel.add(lblMethod, "flowx,cell 0 8 2 1,alignx right");
		
		rbtnDollo = new JRadioButton("Dollo");
		rbtnDollo.setSelected(true);
		rbtnDollo.setHorizontalAlignment(SwingConstants.LEFT);
		rbtnDollo.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnDollo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				MethodToggle(true);
			}
		});
		panel.add(rbtnDollo, "cell 2 8");
		
		rbtnPoly = new JRadioButton("Polymorphism");
		rbtnPoly.setSelected(false);
		rbtnPoly.setHorizontalAlignment(SwingConstants.LEFT);
		rbtnPoly.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnPoly.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				MethodToggle(false);
			}
		});
		panel.add(rbtnPoly, "cell 2 8");

		lblRandOrder = new JLabel("Randomize input order of sequences:");
		lblRandOrder.setFont(new Font("Arial", Font.BOLD, 13));
		lblRandOrder.setHorizontalAlignment(SwingConstants.RIGHT);
		panel.add(lblRandOrder, "flowx,cell 0 9 2 1,alignx right");

		rbtnRandOrderYes = new JRadioButton("Yes");
		rbtnRandOrderYes.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnRandOrderYes.setHorizontalAlignment(SwingConstants.LEFT);
		rbtnRandOrderYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				RandToggle(false);
			}
		});
		rbtnRandOrderYes.setSelected(false);
		panel.add(rbtnRandOrderYes, "cell 2 9");

		rbtnRandOrderNo = new JRadioButton("No, use input order");
		rbtnRandOrderNo.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnRandOrderNo.setHorizontalAlignment(SwingConstants.LEFT);
		rbtnRandOrderNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				RandToggle(true);
			}
		});
		rbtnRandOrderNo.setSelected(true);
		panel.add(rbtnRandOrderNo, "cell 2 9");

		lblRandomSeed = new JLabel("Random number seed:");
		lblRandomSeed.setFont(new Font("Arial", Font.BOLD, 13));
		lblRandomSeed.setHorizontalAlignment(SwingConstants.RIGHT);
		lblRandomSeed.setEnabled(false);
		panel.add(lblRandomSeed, "flowx,cell 0 10 2 1,alignx right");

		txtRandomSeed = new JTextField();
		txtRandomSeed.setText("1");
		txtRandomSeed.setFont(new Font("Arial", Font.PLAIN, 13));
		txtRandomSeed.setEnabled(false);
		txtRandomSeed.setColumns(6);
		panel.add(txtRandomSeed, "cell 2 10");
		
		lblRandOdd = new JLabel("(must be odd)");
		lblRandOdd.setEnabled(false);
		lblRandOdd.setHorizontalAlignment(SwingConstants.LEFT);
		lblRandOdd.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblRandOdd, "cell 2 10");

		lblNumberJumble = new JLabel("Number of times to jumble:");
		lblNumberJumble.setFont(new Font("Arial", Font.BOLD, 13));
		lblNumberJumble.setHorizontalAlignment(SwingConstants.RIGHT);
		lblNumberJumble.setEnabled(false);
		panel.add(lblNumberJumble, "flowx,cell 0 11 2 1,alignx right");

		txtNumberJumble = new JTextField();
		txtNumberJumble.setText("1");
		txtNumberJumble.setFont(new Font("Arial", Font.PLAIN, 13));
		txtNumberJumble.setEnabled(false);
		txtNumberJumble.setColumns(6);
		panel.add(txtNumberJumble, "cell 2 11");

		lblUseThreshold = new JLabel("Use Threshold parsimony:");
		lblUseThreshold.setFont(new Font("Arial", Font.BOLD, 13));
		lblUseThreshold.setHorizontalAlignment(SwingConstants.RIGHT);
		panel.add(lblUseThreshold, "flowx,cell 0 12 2 1,alignx right");

		rbtnThresholdYes = new JRadioButton("Yes");
		rbtnThresholdYes.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnThresholdYes.setHorizontalAlignment(SwingConstants.LEFT);
		rbtnThresholdYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				ThresholdToggle(true);
			}
		});
		rbtnThresholdYes.setSelected(false);
		panel.add(rbtnThresholdYes, "cell 2 12");

		rbtnThresholdNo = new JRadioButton("No, use ordinary parsimony");
		rbtnThresholdNo.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnThresholdNo.setHorizontalAlignment(SwingConstants.LEFT);
		rbtnThresholdNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				ThresholdToggle(false);
			}
		});
		rbtnThresholdNo.setSelected(true);
		panel.add(rbtnThresholdNo, "cell 2 12");

		lblThresholdValue = new JLabel("Threshold value:");
		lblThresholdValue.setFont(new Font("Arial", Font.BOLD, 13));
		lblThresholdValue.setHorizontalAlignment(SwingConstants.RIGHT);
		lblThresholdValue.setEnabled(false);
		panel.add(lblThresholdValue, "flowx,cell 0 13 2 1,alignx right");

		txtThresholdValue = new JTextField();
		txtThresholdValue.setText("1.0");
		txtThresholdValue.setFont(new Font("Arial", Font.PLAIN, 13));
		txtThresholdValue.setEnabled(false);
		txtThresholdValue.setColumns(6);
		panel.add(txtThresholdValue, "cell 2 13");
		
		lblAncStates = new JLabel("Read ancestral states:");
		lblAncStates.setHorizontalAlignment(SwingConstants.RIGHT);
		lblAncStates.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblAncStates, "flowx,cell 0 14 2 1,alignx right");
		
		rbtnAncStatesYes = new JRadioButton("Yes");
		rbtnAncStatesYes.setSelected(false);
		rbtnAncStatesYes.setHorizontalAlignment(SwingConstants.LEFT);
		rbtnAncStatesYes.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnAncStatesYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				ReadAncestorsToggle(true);
			}
		});
		panel.add(rbtnAncStatesYes, "cell 2 14");
		
		rbtnAncStatesNo = new JRadioButton("No");
		rbtnAncStatesNo.setSelected(true);
		rbtnAncStatesNo.setHorizontalAlignment(SwingConstants.LEFT);
		rbtnAncStatesNo.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnAncStatesNo.setBounds(343, 373, 73, 23);
		rbtnAncStatesNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				ReadAncestorsToggle(false);
			}
		});
		panel.add(rbtnAncStatesNo, "cell 2 14");

		lblSitesWeight = new JLabel("Sites weighted:");
		lblSitesWeight.setFont(new Font("Arial", Font.BOLD, 13));
		lblSitesWeight.setHorizontalAlignment(SwingConstants.RIGHT);
		panel.add(lblSitesWeight, "flowx,cell 0 15 2 1,alignx right");

		rbtnSitesYes = new JRadioButton("Yes");
		rbtnSitesYes.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnSitesYes.setHorizontalAlignment(SwingConstants.LEFT);
		rbtnSitesYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				SiteToggle(true);
			}
		});
		rbtnSitesYes.setSelected(false);
		panel.add(rbtnSitesYes, "cell 2 15");

		rbtnSitesNo = new JRadioButton("No");
		rbtnSitesNo.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnSitesNo.setHorizontalAlignment(SwingConstants.LEFT);
		rbtnSitesNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				SiteToggle(false);
			}
		});
		rbtnSitesNo.setSelected(true);
		panel.add(rbtnSitesNo, "cell 2 15");

		lblMultData = new JLabel("Analyze multiple data sets:");
		lblMultData.setFont(new Font("Arial", Font.BOLD, 13));
		lblMultData.setHorizontalAlignment(SwingConstants.RIGHT);
		panel.add(lblMultData, "flowx,cell 0 16 2 1,alignx right");

		rbtnMultDataYes = new JRadioButton("Yes");
		rbtnMultDataYes.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnMultDataYes.setHorizontalAlignment(SwingConstants.LEFT);
		rbtnMultDataYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				MultToggle(true);
			}
		});
		rbtnMultDataYes.setSelected(false);
		panel.add(rbtnMultDataYes, "cell 2 16");

		rbtnMultDataNo = new JRadioButton("No");
		rbtnMultDataNo.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnMultDataNo.setHorizontalAlignment(SwingConstants.LEFT);
		rbtnMultDataNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				MultToggle(false);
			}
		});
		rbtnMultDataNo.setSelected(true);
		panel.add(rbtnMultDataNo, "cell 2 16");

		lblMultDataWeight = new JLabel(
				"Multiple data sets or multiple weights:");
		lblMultDataWeight.setFont(new Font("Arial", Font.BOLD, 13));
		lblMultDataWeight.setHorizontalAlignment(SwingConstants.RIGHT);
		lblMultDataWeight.setEnabled(false);
		panel.add(lblMultDataWeight, "flowx,cell 0 17 2 1,alignx right");

		rbtnDataSet = new JRadioButton("Data sets");
		rbtnDataSet.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnDataSet.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				DataWeightToggle(true);
			}
		});
		rbtnDataSet.setSelected(true);
		rbtnDataSet.setEnabled(false);
		panel.add(rbtnDataSet, "cell 2 17");

		rbtnWeights = new JRadioButton("Weights");
		rbtnWeights.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnWeights.setHorizontalAlignment(SwingConstants.LEFT);
		rbtnWeights.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				DataWeightToggle(false);
			}
		});
		rbtnWeights.setSelected(false);
		rbtnWeights.setEnabled(false);
		panel.add(rbtnWeights, "cell 2 17");
		
		lblNSets = new JLabel("Number:\r\n");
		lblNSets.setHorizontalAlignment(SwingConstants.RIGHT);
		lblNSets.setFont(new Font("Arial", Font.BOLD, 13));
		lblNSets.setEnabled(false);
		panel.add(lblNSets, "cell 2 17");
		
		txtNSets = new JTextField();
		txtNSets.setText("1");
		txtNSets.setFont(new Font("Arial", Font.PLAIN, 13));
		txtNSets.setEnabled(false);
		txtNSets.setColumns(6);
		panel.add(txtNSets, "cell 2 17");

		lblPrintData = new JLabel("Print out the data at start of run:");
		lblPrintData.setFont(new Font("Arial", Font.BOLD, 13));
		lblPrintData.setHorizontalAlignment(SwingConstants.RIGHT);
		panel.add(lblPrintData, "flowx,cell 0 18 2 1,alignx right");

		rbtnPrintDataYes = new JRadioButton("Yes");
		rbtnPrintDataYes.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnPrintDataYes.setHorizontalAlignment(SwingConstants.LEFT);
		rbtnPrintDataYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				PrintDataToggle(true);
			}
		});
		rbtnPrintDataYes.setSelected(false);
		panel.add(rbtnPrintDataYes, "cell 2 18");

		rbtnPrintDataNo = new JRadioButton("No");
		rbtnPrintDataNo.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnPrintDataNo.setHorizontalAlignment(SwingConstants.LEFT);
		rbtnPrintDataNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				PrintDataToggle(false);
			}
		});
		rbtnPrintDataNo.setSelected(true);
		panel.add(rbtnPrintDataNo, "cell 2 18");

		lblPrintTree = new JLabel("Print out tree:");
		lblPrintTree.setFont(new Font("Arial", Font.BOLD, 13));
		lblPrintTree.setHorizontalAlignment(SwingConstants.RIGHT);
		panel.add(lblPrintTree, "flowx,cell 0 19 2 1,alignx right");

		rbtnPrintTreeYes = new JRadioButton("Yes");
		rbtnPrintTreeYes.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnPrintTreeYes.setHorizontalAlignment(SwingConstants.LEFT);
		rbtnPrintTreeYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				PrintTreeToggle(true);
			}
		});
		rbtnPrintTreeYes.setSelected(true);
		panel.add(rbtnPrintTreeYes, "cell 2 19");

		rbtnPrintTreeNo = new JRadioButton("No");
		rbtnPrintTreeNo.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnPrintTreeNo.setHorizontalAlignment(SwingConstants.LEFT);
		rbtnPrintTreeNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				PrintTreeToggle(false);
			}
		});
		rbtnPrintTreeNo.setSelected(false);
		panel.add(rbtnPrintTreeNo, "cell 2 19");

		lblPrintSteps = new JLabel("Print out steps in each character:");
		lblPrintSteps.setFont(new Font("Arial", Font.BOLD, 13));
		lblPrintSteps.setHorizontalAlignment(SwingConstants.RIGHT);
		panel.add(lblPrintSteps, "flowx,cell 0 20 2 1,alignx right");

		rbtnPrintStepYes = new JRadioButton("Yes");
		rbtnPrintStepYes.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnPrintStepYes.setHorizontalAlignment(SwingConstants.LEFT);
		rbtnPrintStepYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				PrintStepToggle(true);
			}
		});
		rbtnPrintStepYes.setSelected(false);
		panel.add(rbtnPrintStepYes, "cell 2 20");

		rbtnPrintStepNo = new JRadioButton("No");
		rbtnPrintStepNo.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnPrintStepNo.setHorizontalAlignment(SwingConstants.LEFT);
		rbtnPrintStepNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				PrintStepToggle(false);
			}
		});
		rbtnPrintStepNo.setSelected(true);
		panel.add(rbtnPrintStepNo, "cell 2 20");

		lblPrintChar = new JLabel("Print states at all nodes of tree:");
		lblPrintChar.setFont(new Font("Arial", Font.BOLD, 13));
		lblPrintChar.setHorizontalAlignment(SwingConstants.RIGHT);
		panel.add(lblPrintChar, "flowx,cell 0 21 2 1,alignx right");

		rbtnPrintCharYes = new JRadioButton("Yes");
		rbtnPrintCharYes.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnPrintCharYes.setHorizontalAlignment(SwingConstants.LEFT);
		rbtnPrintCharYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				PrintSeqToggle(true);
			}
		});
		rbtnPrintCharYes.setSelected(false);
		panel.add(rbtnPrintCharYes, "cell 2 21");

		rbtnPrintCharNo = new JRadioButton("No");
		rbtnPrintCharNo.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnPrintCharNo.setHorizontalAlignment(SwingConstants.LEFT);
		rbtnPrintCharNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				PrintSeqToggle(false);
			}
		});
		rbtnPrintCharNo.setSelected(true);
		panel.add(rbtnPrintCharNo, "cell 2 21");

		lblWriteTree = new JLabel("Write out trees onto tree file:");
		lblWriteTree.setFont(new Font("Arial", Font.BOLD, 13));
		lblWriteTree.setHorizontalAlignment(SwingConstants.RIGHT);
		panel.add(lblWriteTree, "flowx,cell 0 22 2 1,alignx right");

		rbtnWriteTreeYes = new JRadioButton("Yes");
		rbtnWriteTreeYes.setSelected(true);
		rbtnWriteTreeYes.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnWriteTreeYes.setHorizontalAlignment(SwingConstants.LEFT);
		rbtnWriteTreeYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				WriteTreeToggle(true);
			}
		});
		panel.add(rbtnWriteTreeYes, "cell 2 22");

		rbtnWriteTreeNo = new JRadioButton("No");
		rbtnWriteTreeNo.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnWriteTreeNo.setHorizontalAlignment(SwingConstants.LEFT);
		rbtnWriteTreeNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				WriteTreeToggle(false);
			}
		});
		panel.add(rbtnWriteTreeNo, "cell 2 22");

		lblPrintInd = new JLabel("Display progress:");
		lblPrintInd.setFont(new Font("Arial", Font.BOLD, 13));
		lblPrintInd.setHorizontalAlignment(SwingConstants.RIGHT);
		panel.add(lblPrintInd, "flowx,cell 0 23 2 1,alignx right");

		rbtnPrintIndYes = new JRadioButton("Yes");
		rbtnPrintIndYes.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnPrintIndYes.setHorizontalAlignment(SwingConstants.LEFT);
		rbtnPrintIndYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				PrintIndToggle(true);
			}
		});
		rbtnPrintIndYes.setSelected(true);
		panel.add(rbtnPrintIndYes, "cell 2 23");

		rbtnPrintIndNo = new JRadioButton("No");
		rbtnPrintIndNo.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnPrintIndNo.setHorizontalAlignment(SwingConstants.LEFT);
		rbtnPrintIndNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				PrintIndToggle(false);
			}
		});
		rbtnPrintIndNo.setSelected(false);
		panel.add(rbtnPrintIndNo, "cell 2 23");

		btnExecute = new JButton("Execute");
		btnExecute.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				inputvals = new DollopData();
				inputvals.infile = txtInputFile.getText();
				inputvals.intree = txtInTree.getText();
				inputvals.weightfile = txtWeightFile.getText();
				inputvals.ancfile = txtAncFile.getText();
				inputvals.outfile = txtOutputFile.getText();
				inputvals.outfileopt = "w";
				inputvals.outtree = txtOutputTree.getText();
				inputvals.outtreeopt = "w";
				inputvals.BestTree = rbtnSearchYes.isSelected();
				inputvals.UseDollo = rbtnDollo.isSelected();
				inputvals.TreeMax = Integer.parseInt(txtNumberOfTrees.getText());
				inputvals.InputOrder = rbtnRandOrderYes.isSelected();
				inputvals.RandNum = Integer.parseInt(txtRandomSeed.getText());
				inputvals.NumJumble = Integer.parseInt(txtNumberJumble.getText());
				inputvals.ThreshDollop = rbtnThresholdYes.isSelected();
				inputvals.ThreshVal = Double.parseDouble(txtThresholdValue.getText());
				inputvals.UseAncStates = rbtnAncStatesYes.isSelected();
				inputvals.SitesWeighted = rbtnSitesYes.isSelected();
				inputvals.AnalyzeMult = rbtnMultDataYes.isSelected();
				inputvals.MultDataset = rbtnDataSet.isSelected();
				inputvals.NumMult = Integer.parseInt(txtNSets.getText());
				inputvals.PrintData = rbtnPrintDataYes.isSelected();
				inputvals.PrintInd = rbtnPrintIndYes.isSelected();
				inputvals.PrintTree = rbtnPrintTreeYes.isSelected();
				inputvals.PrintSteps = rbtnPrintStepYes.isSelected();
				inputvals.PrintSeq = rbtnPrintCharYes.isSelected();
				inputvals.WriteTree = rbtnWriteTreeYes.isSelected();
	
				btnExecute.setEnabled(false);	
				String title = "Dollop Progress";
				String curdir = System.getProperty("user.dir");
				curdir += "/progress.txt";
				File fl = new File(curdir);
				fl.delete();
				inTitle = title;
				inCurdir = curdir;
	
				if (checkInputVals())
				{
					
			  	    Thread dollopThread = new Thread() {
						public void run() {
							runDollopThreads();
						}
			  	    };
			  	  dollopThread.start();
				}
				btnExecute.setEnabled(true);
			}
		});
		btnExecute.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(btnExecute, "cell 2 24,alignx center");

		btnQuit = new JButton("Quit");
		btnQuit.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				if(phylipCall)
				{
					frmDollopControls.dispose();
					return;
				}
				else
				{
					System.exit(0);
				}
			}
		});
		btnQuit.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(btnQuit, "cell 2 24");
	}

	public boolean checkInputVals(){
		
		// check files
		TestFileNames test = new TestFileNames();
		
		if (!test.DuplicateFileNames(inputvals.infile, "Input", inputvals.outfile, "Output"))
		{			
			return false;		
		}

		if((inputvals.WriteTree) && (inputvals.BestTree))
		{
			if (!test.DuplicateFileNames(inputvals.intree, "Input Tree", inputvals.outtree, "Output Tree"))
			{			
				return false;		
			}
			
		}

		if (!test.FileAvailable(inputvals.infile, "Input"))
		{
			return false;
		}
		
		if (inputvals.UseAncStates)
		{
			if (!test.FileAvailable(inputvals.ancfile, "Ancestor States"))
			{
				return false;
			}
		}
		
		if (inputvals.SitesWeighted)
		{
			if (!test.FileAvailable(inputvals.weightfile, "Weights"))
			{
				return false;
			}
		}
		
		if(!inputvals.BestTree)
		{
			if (!test.FileAvailable(inputvals.intree, "Input Tree"))
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

		// check data
		String msg;
		if (inputvals.RandNum % 2 == 0) 
		{
			msg = "Invalid input. Number seed value must be odd.";
			JOptionPane.showMessageDialog(null, msg, "Error", JOptionPane.ERROR_MESSAGE);
			return false;
		}
		if (inputvals.ThreshVal < 1) 
		{
			msg = "Dollopimony threshold value must be >= 1.0";
			JOptionPane.showMessageDialog(null, msg, "Error", JOptionPane.ERROR_MESSAGE);
			return false;
		}
		return true;
	}
	
	protected void runDollopThreads() {
    	try
    	{
    		// see if library exists
    		Native.loadLibrary("dollop", Dollop.class);
		}
    	catch(UnsatisfiedLinkError e)
    	{
     		new TestFileNames().LibraryMissing("Dollop");
    		return;
    	}
		try 
		{
	  	    Thread dollopRunThread = new Thread() {
		  	      public void run() {
		  	    	  
		  			// at this point we hook into the C code
		  			Dollop dollop = (Dollop) Native.loadLibrary("dollop", Dollop.class);
		  			dollop.dollop(
		  				inputvals.infile,
		  				inputvals.intree,
		  				inputvals.weightfile,
		  				inputvals.ancfile,
		  				inputvals.outfile,
		  				inputvals.outfileopt,
		  				inputvals.outtree,
		  				inputvals.outtreeopt,
		  				inputvals.BestTree,
		  				inputvals.UseDollo,
		  				inputvals.TreeMax,
		  				inputvals.InputOrder,
		  				inputvals.RandNum,
		  				inputvals.NumJumble,
		  				inputvals.ThreshDollop,
		  				inputvals.ThreshVal,
		  				inputvals.UseAncStates,
		  				inputvals.SitesWeighted,
		  				inputvals.AnalyzeMult,
		  				inputvals.MultDataset,
		  				inputvals.NumMult,
		  				inputvals.PrintData,
		  				inputvals.PrintInd,
		  				inputvals.PrintTree,
		  				inputvals.PrintSteps,
		  				inputvals.PrintSeq,
		  				inputvals.WriteTree);
			  	    };
	  	    };
	  	    dollopRunThread.start();

	  	    if (inputvals.PrintInd)
	  	    {
		  	  	do
		  	  	{
					Thread.sleep(1000);
					updateProgress();
				} while (dollopRunThread.isAlive());
	  	    }
		} 
		catch (InterruptedException e) {
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
