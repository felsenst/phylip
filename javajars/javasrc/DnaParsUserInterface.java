package phylip;

import java.awt.EventQueue;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JRadioButton;
import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;
import java.io.File;
import javax.swing.JComboBox;
import javax.swing.JFileChooser;
import javax.swing.JOptionPane;
import javax.swing.JTextField;
import javax.swing.JButton;
import javax.swing.SwingConstants;
import javax.swing.SwingUtilities;

import java.awt.Font;
import java.awt.Color;
import javax.swing.DefaultComboBoxModel;

import utilities.DisplayProgress;
import utilities.TestFileNames;

import com.sun.jna.Library;
import com.sun.jna.Native;

import net.miginfocom.swing.MigLayout;
import javax.swing.JScrollPane;
import java.awt.Dimension;
import javax.swing.JPanel;

public class DnaParsUserInterface{
   public interface DnaPars extends Library {
        public void dnapars(
    		String infile,
    		String intree,
    		String outfile,
    		String outfileopt,
    		String weightfile,
    		String outtree,
    		String outtreeopt,
    		boolean searchbest,
    		String searchopts,
    		int TreeSave,
    		boolean inputorder,
    		int RandNum,
    		int NumJumble,
    		boolean OutRoot,
    		int OutNum,
    		boolean ThreshPars,
    		double ThreshVal,
    		boolean TransPars,
    		boolean SitesWeighted,
    		boolean AnalyzeMult,
    		boolean multdataset,
    		int NumSites,
    		boolean inputseq,
    		boolean PrintData,
    		boolean DotDiff,
    		boolean PrintInd,
    		boolean PrintTree,
    		boolean PrintSteps,
    		boolean PrintSeq,
    		boolean WriteTree);
    }

	public class DnaParsData {
		String infile;
		String intree;
		String outfile;
		String outfileopt;
		String weightfile;
		String outtree;
		String outtreeopt;
		boolean searchbest;
		String searchopts;
		int TreeSave;
		boolean inputorder;
		int RandNum;
		int NumJumble;
		boolean OutRoot;
		int OutNum;
		boolean ThreshPars;
		double ThreshVal;
		boolean TransPars;
		boolean SitesWeighted;
		boolean AnalyzeMult;
		int NumMult;
		boolean multdataset;
		boolean inputseq;
		boolean PrintData;
		boolean DotDiff;
		boolean PrintInd;
		boolean PrintTree;
		boolean PrintSteps;
		boolean PrintSeq;
		boolean WriteTree;
	}

	private DnaParsData inputvals;
	private String inTitle;
	private String inCurdir;
	private DisplayProgress dp;
	private String filedir;
	private boolean ExplicitWgts;
	private boolean phylipCall;

	private JFrame frmDnaParsControls;
	private JLabel lblSearchTree;
	private JRadioButton rdbtnSearchYes;
	private JRadioButton rdbtnSearchNo;
	private JLabel lblSearchOpt;
	private JComboBox cmbxSearchOpts;
	private JLabel lblNumberOfTrees;
	private JTextField txtNumberOfTrees;
	private JLabel lblRandOrder;
	private JRadioButton rdbtnRandOrderYes;
	private JRadioButton rdbtnRandOrderNo;
	private JLabel lblRandomSeed;
	private JTextField txtRandomSeed;
	private JLabel lblNumberJumble;
	private JTextField txtNumberJumble;
	private JLabel lblOutRoot;
	private JRadioButton rdbtnOutYes;
	private JRadioButton rdbtnOutNo;
	private JTextField txtOutNumber;
	private JLabel lblOutNumber;
	private JLabel lblUseThreshold;
	private JRadioButton rdbtnThresholdYes;
	private JRadioButton rdbtnThresholdNo;
	private JLabel lblThresholdValue;
	private JTextField txtThresholdValue;
	private JLabel lblUseTransPars;
	private JRadioButton rdbtnUseTransParsYes;
	private JRadioButton rdbtnUseTransParsNo;
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
	private JLabel lblPrintTree;
	private JRadioButton rdbtnPrintTreeYes;
	private JRadioButton rdbtnPrintTreeNo;
	private JLabel lblPrintSteps;
	private JLabel lblDotDiff;
	private JRadioButton rdbtnDotDiffYes;
	private JRadioButton rdbtnDotDiffNo;
	private JRadioButton rdbtnPrintStepYes;
	private JRadioButton rdbtnPrintStepNo;
	private JLabel lblPrintSeq;
	private JRadioButton rdbtnPrintSeqYes;
	private JRadioButton rdbtnPrintSeqNo;
	private JLabel lblWriteTree;
	private JRadioButton rdbtnWriteTreeYes;
	private JRadioButton rdbtnWriteTreeNo;
	private JTextField txtInputFile;
	private JButton btnInputFile;
	private JButton btnExecute;
	private JButton btnQuit;
	private JTextField txtInputTree;
	private JTextField txtOutputFile;
	private JButton btnInputTree;
	private JButton btnOutputFile;
	private JButton btnWeightFile;
	private JTextField txtWeightFile;
	private JButton btnOutputTree;
	private JTextField txtOutputTree;
	private JLabel lblNSets;
	private JTextField txtNSets;
	private JLabel lblRandOdd;
	
	private JScrollPane scrollPane;
	private JPanel panel;

	/**
	 * Launch the application.
	 */
	public static void main(final String[] args) {
		EventQueue.invokeLater(new Runnable() {
			public void run() {
				try {
					DnaParsUserInterface window = new DnaParsUserInterface(args);
					window.frmDnaParsControls.setVisible(true);
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

		int option = fileChooser.showOpenDialog(frmDnaParsControls.getRootPane());
		if (option == JFileChooser.APPROVE_OPTION) {
			File selectedFile = fileChooser.getSelectedFile();
			filedir = fileChooser.getCurrentDirectory().getAbsolutePath();
			file.setText(selectedFile.getPath());
		}
	}

	protected void SearchToggle(boolean isSearch) {
		if (isSearch) {
			rdbtnSearchYes.setSelected(true);
			rdbtnSearchNo.setSelected(false);
			lblSearchOpt.setEnabled(true);
			cmbxSearchOpts.setEnabled(true);
			lblNumberOfTrees.setEnabled(true);
			txtNumberOfTrees.setEnabled(true);
			lblRandOrder.setEnabled(true);
			rdbtnRandOrderYes.setEnabled(true);
			rdbtnRandOrderNo.setEnabled(true);
			txtRandomSeed.setEditable(true);
			txtNumberJumble.setEditable(true);
			lblOutRoot.setEnabled(true);
			rdbtnOutYes.setEnabled(true);
			rdbtnOutNo.setEnabled(true);
			btnInputTree.setEnabled(false);
			txtInputTree.setEnabled(false);
		} else {
			rdbtnSearchYes.setSelected(false);
			rdbtnSearchNo.setSelected(true);
			lblSearchOpt.setEnabled(false);
			cmbxSearchOpts.setEnabled(false);
			lblNumberOfTrees.setEnabled(false);
			txtNumberOfTrees.setEnabled(false);
			lblRandOrder.setEnabled(false);
			rdbtnRandOrderYes.setEnabled(false);
			rdbtnRandOrderNo.setEnabled(false);
			lblRandomSeed.setEnabled(false);
			txtRandomSeed.setEditable(false);
			lblNumberJumble.setEnabled(false);
			txtNumberJumble.setEditable(false);
			lblOutRoot.setEnabled(false);
			rdbtnOutYes.setEnabled(false);
			rdbtnOutNo.setEnabled(false);
			btnInputTree.setEnabled(true);
			txtInputTree.setEnabled(true);
		}

	}
	
	protected void RandToggle(boolean isRand) {
		if (isRand) {
			rdbtnRandOrderYes.setSelected(false);
			rdbtnRandOrderNo.setSelected(true);
			lblRandomSeed.setEnabled(false);
			txtRandomSeed.setEnabled(false);
			lblNumberJumble.setEnabled(false);
			txtNumberJumble.setEnabled(false);
			lblRandOdd.setEnabled(false);
		} else {
			rdbtnRandOrderYes.setSelected(true);
			rdbtnRandOrderNo.setSelected(false);
			lblRandomSeed.setEnabled(true);
			txtRandomSeed.setEnabled(true);
			lblNumberJumble.setEnabled(true);
			txtNumberJumble.setEnabled(true);
			lblRandOdd.setEnabled(true);
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

	protected void ThresholdToggle(boolean isThresh) {
		if (isThresh) {
			rdbtnThresholdYes.setSelected(true);
			rdbtnThresholdNo.setSelected(false);
			lblThresholdValue.setEnabled(true);
			txtThresholdValue.setEnabled(true);
		} else {
			rdbtnThresholdYes.setSelected(false);
			rdbtnThresholdNo.setSelected(true);
			lblThresholdValue.setEnabled(false);
			txtThresholdValue.setEnabled(false);
		}
	}

	protected void TransToggle(boolean isyes) {
		if (isyes) {
			rdbtnUseTransParsYes.setSelected(true);
			rdbtnUseTransParsNo.setSelected(false);
		} else {
			rdbtnUseTransParsYes.setSelected(false);
			rdbtnUseTransParsNo.setSelected(true);
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

	protected void PrintTreeToggle(boolean isTree) {
		if (isTree) {
			rdbtnPrintTreeYes.setSelected(true);
			rdbtnPrintTreeNo.setSelected(false);
		} else {
			rdbtnPrintTreeYes.setSelected(false);
			rdbtnPrintTreeNo.setSelected(true);
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

	protected void PrintStepToggle(boolean isStep) {
		if (isStep) {
			rdbtnPrintStepYes.setSelected(true);
			rdbtnPrintStepNo.setSelected(false);
		} else {
			rdbtnPrintStepYes.setSelected(false);
			rdbtnPrintStepNo.setSelected(true);
		}
	}

	protected void PrintSeqToggle(boolean isSeq) {
		if (isSeq) {
			rdbtnPrintSeqYes.setSelected(true);
			rdbtnPrintSeqNo.setSelected(false);
		} else {
			rdbtnPrintSeqYes.setSelected(false);
			rdbtnPrintSeqNo.setSelected(true);
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
	
	/**
	 * Create the application.
	 */
	public DnaParsUserInterface(String[] args) {
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
		
		frmDnaParsControls = new JFrame();
		frmDnaParsControls.setBackground(new Color(204, 255, 255));
		frmDnaParsControls.setTitle("Dnapars");
		frmDnaParsControls.setBounds(100, 50, 1050, 620);
		frmDnaParsControls.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frmDnaParsControls.setPreferredSize(new Dimension(frmDnaParsControls.getBounds().width, frmDnaParsControls.getBounds().height));
	
		scrollPane = new JScrollPane();
		scrollPane.setPreferredSize(frmDnaParsControls.getPreferredSize());
		frmDnaParsControls.getContentPane().add(scrollPane);
		
		panel = new JPanel();
		panel.setPreferredSize(frmDnaParsControls.getPreferredSize());
		scrollPane.setViewportView(panel);
		panel.setLayout(new MigLayout("", "[pref!,grow][30.00,grow][pref!,grow][pref!,grow][pref!,grow]", "[][][][]"));

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
		panel.add(txtInputFile, "cell 1 0 4 1,growx");
		
		btnInputTree = new JButton("Input Tree");
		btnInputTree.setEnabled(false);
		btnInputTree.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				ChooseFile(txtInputTree);
			}
		});	
		btnInputTree.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(btnInputTree, "cell 0 1,growx");
		
		txtInputTree = new JTextField();
		txtInputTree.setEnabled(false);
		txtInputTree.setFont(new Font("Arial", Font.PLAIN, 13));
		txtInputTree.setText("intree");
		txtInputTree.setBounds(159, 36, 870, 20);
		panel.add(txtInputTree, "cell 1 1 4 1,growx");
		
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
		txtWeightFile.setBounds(159, 61, 870, 20);
		panel.add(txtWeightFile, "cell 1 2 4 1,growx");
		
		btnOutputTree = new JButton("Output Tree");
		btnOutputTree.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				ChooseFile(txtOutputTree);
			}
		});
		btnOutputTree.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(btnOutputTree, "cell 0 4,growx");
		
		txtOutputTree = new JTextField();
		txtOutputTree.setText("outtree");
		txtOutputTree.setFont(new Font("Arial", Font.PLAIN, 13));
		txtOutputTree.setBounds(159, 111, 870, 20);
		panel.add(txtOutputTree, "cell 1 4 4 1,growx");
	
		btnOutputFile = new JButton("Output File");
		btnOutputFile.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				ChooseFile(txtOutputFile);
			}
		});
		btnOutputFile.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(btnOutputFile, "cell 0 3,growx");
		
		txtOutputFile = new JTextField();
		txtOutputFile.setFont(new Font("Arial", Font.PLAIN, 13));
		txtOutputFile.setText("outfile");
		txtOutputFile.setBounds(159, 86, 870, 20);
		panel.add(txtOutputFile, "cell 1 3 4 1,growx");

		lblSearchTree = new JLabel("Search for best tree:");
		lblSearchTree.setFont(new Font("Arial", Font.BOLD, 13));
		lblSearchTree.setHorizontalAlignment(SwingConstants.RIGHT);
		lblSearchTree.setBounds(78, 150, 195, 14);
		panel.add(lblSearchTree, "flowx,cell 0 5 2 1,alignx right");

		rdbtnSearchYes = new JRadioButton("Yes");
		rdbtnSearchYes.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnSearchYes.setHorizontalAlignment(SwingConstants.LEFT);
		rdbtnSearchYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				SearchToggle(true);
			}
		});
		rdbtnSearchYes.setSelected(true);
		panel.add(rdbtnSearchYes, "cell 2 5");

		rdbtnSearchNo = new JRadioButton("No, use user trees in input file");
		rdbtnSearchNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnSearchNo.setHorizontalAlignment(SwingConstants.LEFT);
		rdbtnSearchNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				SearchToggle(false);
			}
		});
		rdbtnSearchNo.setSelected(false);
		panel.add(rdbtnSearchNo, "cell 2 5");

		lblSearchOpt = new JLabel("Method:");
		lblSearchOpt.setFont(new Font("Arial", Font.BOLD, 13));
		lblSearchOpt.setHorizontalAlignment(SwingConstants.RIGHT);
		panel.add(lblSearchOpt, "flowx,cell 0 6 2 1,alignx right");

		cmbxSearchOpts = new JComboBox();
		cmbxSearchOpts.setModel(new DefaultComboBoxModel(new String[] {"More thorough search", "Less thorough search", "Rearrange on just one best tree"}));
		cmbxSearchOpts.setFont(new Font("Arial", Font.PLAIN, 13));
		panel.add(cmbxSearchOpts, "cell 2 6,growx");

		lblNumberOfTrees = new JLabel("Number of trees to save:");
		lblNumberOfTrees.setFont(new Font("Arial", Font.BOLD, 13));
		lblNumberOfTrees.setHorizontalAlignment(SwingConstants.RIGHT);
		panel.add(lblNumberOfTrees, "flowx,cell 0 7 2 1,alignx right");

		txtNumberOfTrees = new JTextField();
		txtNumberOfTrees.setFont(new Font("Arial", Font.PLAIN, 13));
		txtNumberOfTrees.setText("10000");
		txtNumberOfTrees.setColumns(6);
		panel.add(txtNumberOfTrees, "cell 2 7");

		lblRandOrder = new JLabel("Randomize input order of sequences:");
		lblRandOrder.setFont(new Font("Arial", Font.BOLD, 13));
		lblRandOrder.setHorizontalAlignment(SwingConstants.RIGHT);
		panel.add(lblRandOrder, "flowx,cell 0 8 2 1,alignx right");

		rdbtnRandOrderYes = new JRadioButton("Yes");
		rdbtnRandOrderYes.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnRandOrderYes.setHorizontalAlignment(SwingConstants.LEFT);
		rdbtnRandOrderYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				RandToggle(false);
			}
		});
		rdbtnRandOrderYes.setSelected(false);
		panel.add(rdbtnRandOrderYes, "cell 2 8");

		rdbtnRandOrderNo = new JRadioButton("No, use input order");
		rdbtnRandOrderNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnRandOrderNo.setHorizontalAlignment(SwingConstants.LEFT);
		rdbtnRandOrderNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				RandToggle(true);
			}
		});
		rdbtnRandOrderNo.setSelected(true);
		panel.add(rdbtnRandOrderNo, "cell 2 8");

		lblRandomSeed = new JLabel("Random number seed:");
		lblRandomSeed.setFont(new Font("Arial", Font.BOLD, 13));
		lblRandomSeed.setHorizontalAlignment(SwingConstants.RIGHT);
		lblRandomSeed.setEnabled(false);
		panel.add(lblRandomSeed, "flowx,cell 0 9 2 1,alignx right");

		txtRandomSeed = new JTextField();
		txtRandomSeed.setText("1");
		txtRandomSeed.setFont(new Font("Arial", Font.PLAIN, 13));
		txtRandomSeed.setEnabled(false);
		txtRandomSeed.setColumns(6);
		panel.add(txtRandomSeed, "cell 2 9");
		
		lblRandOdd = new JLabel("(must be odd)");
		lblRandOdd.setEnabled(false);
		lblRandOdd.setHorizontalAlignment(SwingConstants.LEFT);
		lblRandOdd.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblRandOdd, "cell 2 9");

		lblNumberJumble = new JLabel("Number of times to jumble:");
		lblNumberJumble.setFont(new Font("Arial", Font.BOLD, 13));
		lblNumberJumble.setHorizontalAlignment(SwingConstants.RIGHT);
		lblNumberJumble.setEnabled(false);
		panel.add(lblNumberJumble, "flowx,cell 0 10 2 1,alignx right");

		txtNumberJumble = new JTextField();
		txtNumberJumble.setText("1");
		txtNumberJumble.setFont(new Font("Arial", Font.PLAIN, 13));
		txtNumberJumble.setEnabled(false);
		txtNumberJumble.setColumns(6);
		panel.add(txtNumberJumble, "cell 2 10");

		lblOutRoot = new JLabel("Outgroup Root:");
		lblOutRoot.setFont(new Font("Arial", Font.BOLD, 13));
		lblOutRoot.setHorizontalAlignment(SwingConstants.RIGHT);
		panel.add(lblOutRoot, "flowx,cell 0 11 2 1,alignx right");

		rdbtnOutYes = new JRadioButton("Yes");
		rdbtnOutYes.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnOutYes.setHorizontalAlignment(SwingConstants.LEFT);
		rdbtnOutYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				OutToggle(true);
			}
		});
		rdbtnOutYes.setSelected(false);
		panel.add(rdbtnOutYes, "cell 2 11");

		rdbtnOutNo = new JRadioButton("No, use as outgroup species");
		rdbtnOutNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnOutNo.setHorizontalAlignment(SwingConstants.LEFT);
		rdbtnOutNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				OutToggle(false);
			}
		});
		rdbtnOutNo.setSelected(true);
		panel.add(rdbtnOutNo, "cell 2 11");

		lblOutNumber = new JLabel("Number of the outgroup:\r\n");
		lblOutNumber.setFont(new Font("Arial", Font.BOLD, 13));
		lblOutNumber.setHorizontalAlignment(SwingConstants.RIGHT);
		lblOutNumber.setEnabled(false);
		panel.add(lblOutNumber, "flowx,cell 0 12 2 1,alignx right");

		txtOutNumber = new JTextField();
		txtOutNumber.setText("1");
		txtOutNumber.setFont(new Font("Arial", Font.PLAIN, 13));
		txtOutNumber.setEnabled(false);
		txtOutNumber.setColumns(6);
		panel.add(txtOutNumber, "cell 2 12");

		lblUseThreshold = new JLabel("Use Threshold parsimony:");
		lblUseThreshold.setFont(new Font("Arial", Font.BOLD, 13));
		lblUseThreshold.setHorizontalAlignment(SwingConstants.RIGHT);
		panel.add(lblUseThreshold, "flowx,cell 0 13 2 1,alignx right");

		rdbtnThresholdYes = new JRadioButton("Yes");
		rdbtnThresholdYes.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnThresholdYes.setHorizontalAlignment(SwingConstants.LEFT);
		rdbtnThresholdYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				ThresholdToggle(true);
			}
		});
		rdbtnThresholdYes.setSelected(false);
		panel.add(rdbtnThresholdYes, "cell 2 13");

		rdbtnThresholdNo = new JRadioButton("No, use ordinary parsimony");
		rdbtnThresholdNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnThresholdNo.setHorizontalAlignment(SwingConstants.LEFT);
		rdbtnThresholdNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				ThresholdToggle(false);
			}
		});
		rdbtnThresholdNo.setSelected(true);
		panel.add(rdbtnThresholdNo, "cell 2 13");

		lblThresholdValue = new JLabel("Threshold value:");
		lblThresholdValue.setFont(new Font("Arial", Font.BOLD, 13));
		lblThresholdValue.setHorizontalAlignment(SwingConstants.RIGHT);
		lblThresholdValue.setEnabled(false);
		panel.add(lblThresholdValue, "flowx,cell 0 14 2 1,alignx right");

		txtThresholdValue = new JTextField();
		txtThresholdValue.setText("1.0");
		txtThresholdValue.setFont(new Font("Arial", Font.PLAIN, 13));
		txtThresholdValue.setEnabled(false);
		txtThresholdValue.setColumns(6);
		panel.add(txtThresholdValue, "cell 2 14");

		lblUseTransPars = new JLabel("Use Transversion parsimony:");
		lblUseTransPars.setFont(new Font("Arial", Font.BOLD, 13));
		lblUseTransPars.setHorizontalAlignment(SwingConstants.RIGHT);
		panel.add(lblUseTransPars, "flowx,cell 0 15 2 1,alignx right");

		rdbtnUseTransParsYes = new JRadioButton("Yes");
		rdbtnUseTransParsYes.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnUseTransParsYes.setHorizontalAlignment(SwingConstants.LEFT);
		rdbtnUseTransParsYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				TransToggle(true);
			}
		});
		rdbtnUseTransParsYes.setSelected(false);
		panel.add(rdbtnUseTransParsYes, "cell 2 15");

		rdbtnUseTransParsNo = new JRadioButton("No, count all steps");
		rdbtnUseTransParsNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnUseTransParsNo.setHorizontalAlignment(SwingConstants.LEFT);
		rdbtnUseTransParsNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				TransToggle(false);
			}
		});
		rdbtnUseTransParsNo.setSelected(true);
		panel.add(rdbtnUseTransParsNo, "cell 2 15");

		// **** column 2 **** //
		lblSitesWeight = new JLabel("Sites weighted:");
		lblSitesWeight.setFont(new Font("Arial", Font.BOLD, 13));
		lblSitesWeight.setHorizontalAlignment(SwingConstants.RIGHT);
		panel.add(lblSitesWeight, "flowx,cell 0 16 2 1,alignx right");

		rdbtnSitesYes = new JRadioButton("Yes");
		rdbtnSitesYes.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnSitesYes.setHorizontalAlignment(SwingConstants.LEFT);
		rdbtnSitesYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				SiteToggle(true);
			}
		});
		rdbtnSitesYes.setSelected(false);
		panel.add(rdbtnSitesYes, "cell 2 16");

		rdbtnSitesNo = new JRadioButton("No");
		rdbtnSitesNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnSitesNo.setHorizontalAlignment(SwingConstants.LEFT);
		rdbtnSitesNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				SiteToggle(false);
			}
		});
		rdbtnSitesNo.setSelected(true);
		panel.add(rdbtnSitesNo, "cell 2 16");

		lblMultData = new JLabel("Analyze multiple data sets:");
		lblMultData.setFont(new Font("Arial", Font.BOLD, 13));
		lblMultData.setHorizontalAlignment(SwingConstants.RIGHT);
		panel.add(lblMultData, "cell 3 5,alignx right");

		rdbtnMultDataYes = new JRadioButton("Yes");
		rdbtnMultDataYes.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnMultDataYes.setHorizontalAlignment(SwingConstants.LEFT);
		rdbtnMultDataYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				MultToggle(true);
			}
		});
		rdbtnMultDataYes.setSelected(false);
		panel.add(rdbtnMultDataYes, "cell 4 5");

		rdbtnMultDataNo = new JRadioButton("No");
		rdbtnMultDataNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnMultDataNo.setHorizontalAlignment(SwingConstants.LEFT);
		rdbtnMultDataNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				MultToggle(false);
			}
		});
		rdbtnMultDataNo.setSelected(true);
		panel.add(rdbtnMultDataNo, "cell 4 5");

		lblMultDataWeight = new JLabel("Multiple data sets or multiple weights:\r\n");
		lblMultDataWeight.setFont(new Font("Arial", Font.BOLD, 13));
		lblMultDataWeight.setHorizontalAlignment(SwingConstants.RIGHT);
		lblMultDataWeight.setEnabled(false);
		panel.add(lblMultDataWeight, "cell 3 6,alignx right");

		rdbtnDataSet = new JRadioButton("Data sets");
		rdbtnDataSet.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnDataSet.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				DataWeightToggle(true);
			}
		});
		rdbtnDataSet.setSelected(true);
		rdbtnDataSet.setEnabled(false);
		panel.add(rdbtnDataSet, "cell 4 6");

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
		panel.add(rdbtnWeights, "cell 4 6");
		
		lblNSets = new JLabel("Number:\r\n");
		lblNSets.setHorizontalAlignment(SwingConstants.RIGHT);
		lblNSets.setFont(new Font("Arial", Font.BOLD, 13));
		lblNSets.setEnabled(false);
		panel.add(lblNSets, "cell 4 7");
		
		txtNSets = new JTextField();
		txtNSets.setText("1");
		txtNSets.setFont(new Font("Arial", Font.PLAIN, 13));
		txtNSets.setEnabled(false);
		txtNSets.setColumns(6);
		panel.add(txtNSets, "cell 4 7");

		lblInputSeq = new JLabel("Input sequences:");
		lblInputSeq.setEnabled(false);
		lblInputSeq.setFont(new Font("Arial", Font.BOLD, 13));
		lblInputSeq.setHorizontalAlignment(SwingConstants.RIGHT);
		panel.add(lblInputSeq, "cell 3 8,alignx right");

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
		panel.add(rdbtnInputSeqYes, "cell 4 8");

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
		panel.add(rdbtnInputSeqNo, "cell 4 8");

		lblPrintData = new JLabel("Print out the data at start of run:");
		lblPrintData.setFont(new Font("Arial", Font.BOLD, 13));
		lblPrintData.setHorizontalAlignment(SwingConstants.RIGHT);
		panel.add(lblPrintData, "cell 3 9,alignx right");

		rdbtnPrintDataYes = new JRadioButton("Yes");
		rdbtnPrintDataYes.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnPrintDataYes.setHorizontalAlignment(SwingConstants.LEFT);
		rdbtnPrintDataYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				PrintDataToggle(true);
			}
		});
		rdbtnPrintDataYes.setSelected(false);
		panel.add(rdbtnPrintDataYes, "cell 4 9");

		rdbtnPrintDataNo = new JRadioButton("No");
		rdbtnPrintDataNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnPrintDataNo.setHorizontalAlignment(SwingConstants.LEFT);
		rdbtnPrintDataNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				PrintDataToggle(false);
			}
		});
		rdbtnPrintDataNo.setSelected(true);
		panel.add(rdbtnPrintDataNo, "cell 4 9");

		lblDotDiff = new JLabel("Use dot-differencing to display them:");
		lblDotDiff.setFont(new Font("Arial", Font.BOLD, 13));
		lblDotDiff.setHorizontalAlignment(SwingConstants.RIGHT);
		panel.add(lblDotDiff, "cell 3 10,alignx right");

		rdbtnDotDiffYes = new JRadioButton("Yes");
		rdbtnDotDiffYes.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnDotDiffYes.setHorizontalAlignment(SwingConstants.LEFT);
		rdbtnDotDiffYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				DotDiffToggle(true);
			}
		});
		rdbtnDotDiffYes.setSelected(true);
		panel.add(rdbtnDotDiffYes, "cell 4 10");

		rdbtnDotDiffNo = new JRadioButton("No");
		rdbtnDotDiffNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnDotDiffNo.setHorizontalAlignment(SwingConstants.LEFT);
		rdbtnDotDiffNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				DotDiffToggle(false);
			}
		});
		panel.add(rdbtnDotDiffNo, "cell 4 10");

		lblPrintTree = new JLabel("Print out tree:");
		lblPrintTree.setFont(new Font("Arial", Font.BOLD, 13));
		lblPrintTree.setHorizontalAlignment(SwingConstants.RIGHT);
		panel.add(lblPrintTree, "cell 3 11,alignx right");

		rdbtnPrintTreeYes = new JRadioButton("Yes");
		rdbtnPrintTreeYes.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnPrintTreeYes.setHorizontalAlignment(SwingConstants.LEFT);
		rdbtnPrintTreeYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				PrintTreeToggle(true);
			}
		});
		rdbtnPrintTreeYes.setSelected(true);
		panel.add(rdbtnPrintTreeYes, "cell 4 11");

		rdbtnPrintTreeNo = new JRadioButton("No");
		rdbtnPrintTreeNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnPrintTreeNo.setHorizontalAlignment(SwingConstants.LEFT);
		rdbtnPrintTreeNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				PrintTreeToggle(false);
			}
		});
		rdbtnPrintTreeNo.setSelected(false);
		panel.add(rdbtnPrintTreeNo, "cell 4 11");

		lblPrintSteps = new JLabel("Print out steps in each site:");
		lblPrintSteps.setFont(new Font("Arial", Font.BOLD, 13));
		lblPrintSteps.setHorizontalAlignment(SwingConstants.RIGHT);
		panel.add(lblPrintSteps, "cell 3 12,alignx right");

		rdbtnPrintStepYes = new JRadioButton("Yes");
		rdbtnPrintStepYes.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnPrintStepYes.setHorizontalAlignment(SwingConstants.LEFT);
		rdbtnPrintStepYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				PrintStepToggle(true);
			}
		});
		rdbtnPrintStepYes.setSelected(false);
		panel.add(rdbtnPrintStepYes, "cell 4 12");

		rdbtnPrintStepNo = new JRadioButton("No");
		rdbtnPrintStepNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnPrintStepNo.setHorizontalAlignment(SwingConstants.LEFT);
		rdbtnPrintStepNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				PrintStepToggle(false);
			}
		});
		rdbtnPrintStepNo.setSelected(true);
		panel.add(rdbtnPrintStepNo, "cell 4 12");

		lblPrintSeq = new JLabel("Print sequences at all nodes of tree:");
		lblPrintSeq.setFont(new Font("Arial", Font.BOLD, 13));
		lblPrintSeq.setHorizontalAlignment(SwingConstants.RIGHT);
		panel.add(lblPrintSeq, "cell 3 13,alignx right");

		rdbtnPrintSeqYes = new JRadioButton("Yes");
		rdbtnPrintSeqYes.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnPrintSeqYes.setHorizontalAlignment(SwingConstants.LEFT);
		rdbtnPrintSeqYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				PrintSeqToggle(true);
			}
		});
		rdbtnPrintSeqYes.setSelected(false);
		panel.add(rdbtnPrintSeqYes, "cell 4 13");

		rdbtnPrintSeqNo = new JRadioButton("No");
		rdbtnPrintSeqNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnPrintSeqNo.setHorizontalAlignment(SwingConstants.LEFT);
		rdbtnPrintSeqNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				PrintSeqToggle(false);
			}
		});
		rdbtnPrintSeqNo.setSelected(true);
		panel.add(rdbtnPrintSeqNo, "cell 4 13");

		lblWriteTree = new JLabel("Write out trees onto tree file:");
		lblWriteTree.setFont(new Font("Arial", Font.BOLD, 13));
		lblWriteTree.setHorizontalAlignment(SwingConstants.RIGHT);
		panel.add(lblWriteTree, "cell 3 14,alignx right");

		rdbtnWriteTreeYes = new JRadioButton("Yes");
		rdbtnWriteTreeYes.setSelected(true);
		rdbtnWriteTreeYes.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnWriteTreeYes.setHorizontalAlignment(SwingConstants.LEFT);
		rdbtnWriteTreeYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				WriteTreeToggle(true);
			}
		});
		panel.add(rdbtnWriteTreeYes, "cell 4 14");

		rdbtnWriteTreeNo = new JRadioButton("No");
		rdbtnWriteTreeNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnWriteTreeNo.setHorizontalAlignment(SwingConstants.LEFT);
		rdbtnWriteTreeNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				WriteTreeToggle(false);
			}
		});
		panel.add(rdbtnWriteTreeNo, "cell 4 14");

		lblPrintInd = new JLabel("Display progress:");
		lblPrintInd.setFont(new Font("Arial", Font.BOLD, 13));
		lblPrintInd.setHorizontalAlignment(SwingConstants.RIGHT);
		panel.add(lblPrintInd, "cell 3 15,alignx right");

		rdbtnPrintIndYes = new JRadioButton("Yes");
		rdbtnPrintIndYes.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnPrintIndYes.setHorizontalAlignment(SwingConstants.LEFT);
		rdbtnPrintIndYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				PrintIndToggle(true);
			}
		});
		rdbtnPrintIndYes.setSelected(true);
		panel.add(rdbtnPrintIndYes, "cell 4 15");

		rdbtnPrintIndNo = new JRadioButton("No");
		rdbtnPrintIndNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnPrintIndNo.setHorizontalAlignment(SwingConstants.LEFT);
		rdbtnPrintIndNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				PrintIndToggle(false);
			}
		});
		rdbtnPrintIndNo.setSelected(false);
		panel.add(rdbtnPrintIndNo, "cell 4 15");

		btnExecute = new JButton("Execute");
		btnExecute.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				inputvals = new DnaParsData();
				inputvals.infile = txtInputFile.getText();
				inputvals.intree = txtInputTree.getText();
				inputvals.outfile = txtOutputFile.getText();
				inputvals.outfileopt = "w";
				inputvals.weightfile = txtWeightFile.getText();
				inputvals.outtree = txtOutputTree.getText();
				inputvals.outtreeopt = "w";
				inputvals.searchbest = rdbtnSearchYes.isSelected();
				inputvals.searchopts = cmbxSearchOpts.getSelectedItem().toString();
				inputvals.TreeSave = Integer.parseInt(txtNumberOfTrees.getText());
				inputvals.inputorder = rdbtnRandOrderYes.isSelected();
				inputvals.RandNum = Integer.parseInt(txtRandomSeed.getText());
				inputvals.NumJumble = Integer.parseInt(txtNumberJumble.getText());
				inputvals.OutRoot = rdbtnOutYes.isSelected();
				inputvals.OutNum = Integer.parseInt(txtOutNumber.getText());
				inputvals.ThreshPars = rdbtnThresholdYes.isSelected();
				inputvals.ThreshVal = Double.parseDouble(txtThresholdValue.getText());
				inputvals.TransPars = rdbtnUseTransParsYes.isSelected();
				inputvals.SitesWeighted = rdbtnSitesYes.isSelected();
				inputvals.AnalyzeMult = rdbtnMultDataYes.isSelected();
				inputvals.multdataset = rdbtnDataSet.isSelected();
				inputvals.NumMult = Integer.parseInt(txtNSets.getText());
				inputvals.inputseq = rdbtnInputSeqYes.isSelected();
				inputvals.PrintData = rdbtnPrintDataYes.isSelected();
				inputvals.DotDiff = rdbtnDotDiffYes.isSelected();
				inputvals.PrintInd = rdbtnPrintIndYes.isSelected();
				inputvals.PrintTree = rdbtnPrintTreeYes.isSelected();
				inputvals.PrintSteps = rdbtnPrintStepYes.isSelected();
				inputvals.PrintSeq = rdbtnPrintSeqYes.isSelected();
				inputvals.WriteTree = rdbtnWriteTreeYes.isSelected();

				
				btnExecute.setEnabled(false);	
				String title = "Dnapars Progress";
				String curdir = System.getProperty("user.dir");
				curdir += "/progress.txt";
				File fl = new File(curdir);
				fl.delete();
				inTitle = title;
				inCurdir = curdir;
	
				if (checkInputVals())
				{
					
			  	    Thread dnaParsThread = new Thread() {
						public void run() {
							runDnaParsThreads();
						}
			  	    };
			  	    dnaParsThread.start();
				}
				btnExecute.setEnabled(true);
			}
		});
		btnExecute.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(btnExecute, "cell 3 17,alignx right");

		btnQuit = new JButton("Quit");
		btnQuit.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				if(phylipCall)
				{
					frmDnaParsControls.dispose();
					return;
				}
				else
				{
					System.exit(0);
				}
			}
		});
		btnQuit.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(btnQuit, "cell 4 17");
	}
	public boolean checkInputVals(){
		
		// check files
		TestFileNames test = new TestFileNames();
		
		if (!test.DuplicateFileNames(inputvals.infile, "Input", inputvals.outfile, "Output"))
		{			
			return false;		
		}

		if((inputvals.WriteTree) && (!inputvals.searchbest))
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
		
		if (inputvals.OutNum > 5) 
		{
			msg = "Bad outgroup number. Must be in range 1-5.";
			JOptionPane.showMessageDialog(null, msg, "Error", JOptionPane.ERROR_MESSAGE);
			return false;
		}
		if (inputvals.OutNum < 1) 
		{
			msg = "Bad outgroup number. Must be in range 1-5.";
			JOptionPane.showMessageDialog(null, msg, "Error", JOptionPane.ERROR_MESSAGE);
			return false;
		}

		if (inputvals.ThreshVal < 1.0) 
		{
			msg = "Threshold value cannot be less than 1.0";
			JOptionPane.showMessageDialog(null, msg, "Error", JOptionPane.ERROR_MESSAGE);
				return false;
		}
		
		return true;
	}
	
	protected void runDnaParsThreads() {
    	try
    	{
    		// see if library exists
    		Native.loadLibrary("dnapars", DnaPars.class);
		}
    	catch(UnsatisfiedLinkError e)
    	{
     		new TestFileNames().LibraryMissing("DnaPars");
    		return;
    	}
		try 
		{
	  	    Thread dnaParsRunThread = new Thread() {
		  	      public void run() {
		  	    	  
		  			// at this point we hook into the C code
				  DnaPars dnapars = (DnaPars) Native.loadLibrary("dnapars", DnaPars.class);
		          dnapars.dnapars(
		        		inputvals.infile,
		        		inputvals.intree,
		        		inputvals.outfile,
		        		inputvals.outfileopt,
		        		inputvals.weightfile,
		        		inputvals.outtree,
		        		inputvals.outtreeopt,
		        		inputvals.searchbest,
		        		inputvals.searchopts,
		        		inputvals.TreeSave,
		        		inputvals.inputorder,
		        		inputvals.RandNum,
		        		inputvals.NumJumble,
		        		inputvals.OutRoot,
		        		inputvals.OutNum,
		        		inputvals.ThreshPars,
		        		inputvals.ThreshVal,
		        		inputvals.TransPars,
		        		inputvals.SitesWeighted,
		        		inputvals.AnalyzeMult,
		        		inputvals.multdataset,
		        		inputvals.NumMult,
		        		inputvals.inputseq,
		        		inputvals.PrintData,
		        		inputvals.DotDiff,
		        		inputvals.PrintInd,
		        		inputvals.PrintTree,
		        		inputvals.PrintSteps,
		        		inputvals.PrintSeq,
		        		inputvals.WriteTree);
			  	    };
	  	    };
	  	    dnaParsRunThread.start();

	  	    if (inputvals.PrintInd)
	  	    {
		  	  	do
		  	  	{
					Thread.sleep(1000);
					updateProgress();
				} while (dnaParsRunThread.isAlive());
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