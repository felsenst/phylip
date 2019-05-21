package phylip;

import java.awt.EventQueue;
import java.io.File;

import javax.swing.JFrame;
import javax.swing.JTextField;

import javax.swing.JOptionPane;
import javax.swing.JLabel;
import javax.swing.JRadioButton;
import javax.swing.SwingConstants;
import javax.swing.SwingUtilities;

import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;
import javax.swing.JButton;
import javax.swing.JFileChooser;

//import com.sun.jna.Library;
//import com.sun.jna.Native;

//import utilities.DisplayProgress;
//import utilities.TestFileNames;

import java.awt.Font;
import java.awt.Color;

//import net.miginfocom.swing.MigLayout;
import javax.swing.JScrollPane;
import java.awt.Dimension;
import javax.swing.JPanel;
/*
public class CliqueUserInterface {
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
	*/
/*
	public class cliqueData{
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
	}
	*/

	//private cliqueData inputvals;
	private String inTitle;
	private String inCurdir;
	private DisplayProgress dp;
	private String filedir;
	private boolean phylipCall;
	
	public enum Wgttype{NONE, DATA, WEIGHT}

	private JFrame frmCliqueControls;
	private JTextField txtInputFile;
	private JTextField txtOutfile;
	
	private JTextField txtOutRoot;
	
	private JRadioButton rdbtnAncN;
	private JRadioButton rdbtnAncY;
	private JTextField txtAncfile;
	
	private JRadioButton rdbtnFacN;
	private JRadioButton rdbtnFacY;
	private JTextField txtFactfile;
	
	private JRadioButton rdbtnWgtN;
	private JRadioButton rdbtnWgtY;
	private JTextField txtWgtsfile;
	private boolean weightstate;
	
	private JRadioButton rdbtnCliq;
	private JLabel lblSize; 
	private JTextField txtCliqueSize;
	
	private JRadioButton rdbtnBrunN;
	private JRadioButton rdbtnBrunY;
	
	private JRadioButton rdbtnCmatN;
	private JRadioButton rdbtnCmatY;
	
	private JRadioButton rdbtnPtreeN;
	private JRadioButton rdbtnPtreeY;
	
	private JRadioButton rdbtnOtreeN;
	private JRadioButton rdbtnOtreeY;
	private JTextField txtOuttree;
	
	private JLabel lblNumSets;
	private JTextField txtSets;
	
	private JLabel lblMultipleDatasets;
	private JRadioButton rdbtnMdsY;
	private JRadioButton rdbtnMdsN;
	private JRadioButton rdbtnMwsY;
	private JRadioButton rdbtnOgrpY;
	private JRadioButton rdbtnOgrpN;
	private JLabel lblOutGroup;
	private JButton btnInputFile;
	private JButton btnOutFile;
	private JButton btnOutTree;
	private JButton btnFactFile;
	private JButton btnAncStates;
	private JButton btnWgtsFile;
	private JButton btnExecute;
	private JButton btnQuit;
	private JLabel lblPrintInd;
	private JRadioButton rdbtnPrintIndYes;
	private JRadioButton rdbtnPrintIndNo;
	private JLabel lblPrintBefore;
	private JLabel lblPrintMatrix;
	private JLabel lblPrintTree;
	private JLabel lblWriteTree;
	private JLabel lblUseAnc;
	private JLabel lblUseFac;
	private JLabel lblSiteWgt;
	private JLabel lblMinCliq;
	private JLabel lblOutRoot;
	
	private JScrollPane scrollPane;
	private JPanel panel;


	/**
	 * Launch the application.
	 */
	public static void main(final String[] args) {
		EventQueue.invokeLater(new Runnable() {
			public void run() {
				try {
					CliqueUserInterface window = new CliqueUserInterface(args);
					window.frmCliqueControls.setVisible(true);
				} catch (Exception e) {
					e.printStackTrace();
				}
				return;
			}
		});
	}

	/**
	 * Create the application.
	 */
	public CliqueUserInterface(String[] args) {
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
			btnAncStates.setEnabled(true);
			txtAncfile.setEnabled(true);
		}	
		else{
			rdbtnAncY.setSelected(false);
			rdbtnAncN.setSelected(true);
			btnAncStates.setEnabled(false);
			txtAncfile.setEnabled(false);
		}		
	}
	
	protected void FacToggle(boolean isyes) {		
		if (isyes){
			rdbtnFacY.setSelected(true);
			rdbtnFacN.setSelected(false);
			btnFactFile.setEnabled(true);
			txtFactfile.setEnabled(true);
		}	
		else{
			rdbtnFacY.setSelected(false);
			rdbtnFacN.setSelected(true);
			btnFactFile.setEnabled(false);
			txtFactfile.setEnabled(false);
		}		
	}
	
	protected void WgtToggle(boolean isyes) {		
		if (isyes){
			if (rdbtnMwsY.isSelected() == true){
				rdbtnWgtY.setSelected(false);
				String msg = "Multiple data sets take precidence over multiple weights. \"Yes\" cannot be selected.";
				JOptionPane.showMessageDialog(null, msg, "Warning", JOptionPane.WARNING_MESSAGE);
				return;
			}
			rdbtnWgtY.setSelected(true);
			rdbtnWgtN.setSelected(false);
			btnWgtsFile.setEnabled(true);
			txtWgtsfile.setEnabled(true);
			weightstate = true;
		}	
		else{
			rdbtnWgtY.setSelected(false);
			rdbtnWgtN.setSelected(true);
			btnWgtsFile.setEnabled(false);
			txtWgtsfile.setEnabled(false);
			weightstate = false;
		}		
	}
	
	protected void BrunToggle(boolean isyes) {		
		if (isyes){
			rdbtnBrunY.setSelected(true);
			rdbtnBrunN.setSelected(false);
		}	
		else{
			rdbtnBrunY.setSelected(false);
			rdbtnBrunN.setSelected(true);
		}		
	}
	
	protected void CmatToggle(boolean isyes) {		
		if (isyes){
			rdbtnCmatY.setSelected(true);
			rdbtnCmatN.setSelected(false);
		}	
		else{
			rdbtnCmatY.setSelected(false);
			rdbtnCmatN.setSelected(true);
		}		
	}
	
	protected void PtreeToggle(boolean isyes) {		
		if (isyes){
			rdbtnPtreeY.setSelected(true);
			rdbtnPtreeN.setSelected(false);
		}	
		else{
			rdbtnPtreeY.setSelected(false);
			rdbtnPtreeN.setSelected(true);
		}		
	}
	
	protected void OtreeToggle(boolean isyes) {		
		if (isyes){
			rdbtnOtreeY.setSelected(true);
			rdbtnOtreeN.setSelected(false);
			btnOutTree.setEnabled(true);
			txtOuttree.setEnabled(true);
		}	
		else{
			rdbtnOtreeY.setSelected(false);
			rdbtnOtreeN.setSelected(true);
			btnOutTree.setEnabled(false);
			txtOuttree.setEnabled(false);
		}		
	}
	
	protected void CliqSizeToggle(boolean isno) {		
		if (isno){
			rdbtnCliq.setSelected(true);
			lblSize.setEnabled(false);
			txtCliqueSize.setEnabled(false);
		}	
		else{
			rdbtnCliq.setSelected(false);
			lblSize.setEnabled(true);
			txtCliqueSize.setEnabled(true);
		}		
	}
	
	protected void OutGroupToggle(boolean isyes) {		
		if (isyes){
			rdbtnOgrpY.setSelected(true);
			rdbtnOgrpN.setSelected(false);
			lblOutGroup.setText("At species number:");
		}	
		else{
			rdbtnOgrpY.setSelected(false);
			rdbtnOgrpN.setSelected(true);
			lblOutGroup.setText("Outgroup species:");
		}		
	}
	
	protected void MultiSetToggle(Wgttype kind) {		
		if (kind == Wgttype.DATA){
			rdbtnMdsY.setSelected(true);
			rdbtnMdsN.setSelected(false);
			rdbtnMwsY.setSelected(false);
			lblNumSets.setEnabled(true);
			txtSets.setEnabled(true);
			if (weightstate == true){
				rdbtnWgtY.setSelected(true);
				rdbtnWgtN.setSelected(false);
				btnWgtsFile.setEnabled(true);
				txtWgtsfile.setEnabled(true);
			}			
			else{
				btnWgtsFile.setEnabled(false);
				txtWgtsfile.setEnabled(false);					
			}
		}	
		else{
			if (kind == Wgttype.NONE){
				rdbtnMdsY.setSelected(false);
				rdbtnMdsN.setSelected(true);
				rdbtnMwsY.setSelected(false);	
				lblNumSets.setEnabled(false);
				txtSets.setEnabled(false);
				if (weightstate == true){
					rdbtnWgtY.setSelected(true);
					rdbtnWgtN.setSelected(false);
					btnWgtsFile.setEnabled(true);
					txtWgtsfile.setEnabled(true);
				}
				else{
					btnWgtsFile.setEnabled(false);
					txtWgtsfile.setEnabled(false);					
				}
				
			}
			else {
				// (kind == Wgttype.WEIGHT)
				rdbtnMdsY.setSelected(false);
				rdbtnMdsN.setSelected(false);
				rdbtnMwsY.setSelected(true);
				lblNumSets.setEnabled(true);
				txtSets.setEnabled(true);
				
				if (rdbtnWgtY.isSelected() == true){
					weightstate = true;
					rdbtnWgtY.setSelected(false);
					rdbtnWgtN.setSelected(true);
					String msg = "Multiple weight sets take precidence over multiple weights. Multiple weights have been shut off.";
					msg += "\nUncheck Multiple weight sets to restore Multiple weights state.";
					JOptionPane.showMessageDialog(null, msg, "Warning", JOptionPane.WARNING_MESSAGE);
				}
				else {
					weightstate = false;
				}
				btnWgtsFile.setEnabled(true);
				txtWgtsfile.setEnabled(true);
			}
		}		
	}
	
	protected void PrintIndToggle(boolean isInd){
		if(isInd){
			rdbtnPrintIndYes.setSelected(true);
			rdbtnPrintIndNo.setSelected(false);
		}
		else{
			rdbtnPrintIndYes.setSelected(false);
			rdbtnPrintIndNo.setSelected(true);
		}
	}

	protected void ChooseFile(JTextField file) {		
		JFileChooser fileChooser = new JFileChooser( filedir);

		int option = fileChooser.showOpenDialog(frmCliqueControls.getRootPane());
		if (option == JFileChooser.APPROVE_OPTION) {
			File selectedFile = fileChooser.getSelectedFile();
			filedir = fileChooser.getCurrentDirectory().getAbsolutePath();
			file.setText(selectedFile.getPath());
		}	
	}

	/**
	 * Initialize the contents of the frame.
	 */
	private void initialize() {
		filedir = System.getProperty("user.dir");
		
		frmCliqueControls = new JFrame();
		frmCliqueControls.setBackground(new Color(204, 255, 255));
		frmCliqueControls.setFont(new Font("Arial", Font.BOLD, 13));
		frmCliqueControls.setTitle("Clique");
		frmCliqueControls.setBounds(100, 100, 650, 625);
		frmCliqueControls.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frmCliqueControls.setPreferredSize(new Dimension(frmCliqueControls.getBounds().width, frmCliqueControls.getBounds().height));
		
		scrollPane = new JScrollPane();
		scrollPane.setPreferredSize(frmCliqueControls.getPreferredSize());
		frmCliqueControls.getContentPane().add(scrollPane);
		
		panel = new JPanel();
		panel.setPreferredSize(frmCliqueControls.getPreferredSize());
		scrollPane.setViewportView(panel);
		scrollPane.setViewportView(panel);
		panel.setLayout(new MigLayout("", "[pref!,grow][30.00,grow][pref!,grow][pref!,grow]", "[][][][][]"));
				
		
		btnInputFile = new JButton("Input File");
		btnInputFile.setFont(new Font("Arial", Font.BOLD, 13));
		btnInputFile.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				ChooseFile(txtInputFile);
			}
		});
		panel.add(btnInputFile, "cell 0 0,growx");

		txtInputFile = new JTextField();
		txtInputFile.setFont(new Font("Arial", Font.PLAIN, 13));
		txtInputFile.setText("infile");
		panel.add(txtInputFile, "cell 1 0 3 1,growx");
		
		btnAncStates = new JButton("Ancestral States File");
		btnAncStates.setFont(new Font("Arial", Font.BOLD, 13));
		btnAncStates.setEnabled(false);
		btnAncStates.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				ChooseFile(txtAncfile);
			}
		});
		panel.add(btnAncStates, "cell 0 1,growx");
		
		txtAncfile = new JTextField();
		txtAncfile.setFont(new Font("Arial", Font.PLAIN, 13));
		txtAncfile.setEnabled(false);
		txtAncfile.setText("ancfile");
		panel.add(txtAncfile, "cell 1 1 3 1,growx");
	
		btnFactFile = new JButton("Factors File");
		btnFactFile.setFont(new Font("Arial", Font.BOLD, 13));
		btnFactFile.setEnabled(false);
		btnFactFile.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				ChooseFile(txtFactfile);
			}
		});
		panel.add(btnFactFile, "cell 0 2,growx");
		
		txtFactfile = new JTextField();
		txtFactfile.setFont(new Font("Arial", Font.PLAIN, 13));
		txtFactfile.setEnabled(false);
		txtFactfile.setText("factfile");
		panel.add(txtFactfile, "cell 1 2 3 1,growx");
		
		btnWgtsFile = new JButton("Weights File");
		btnWgtsFile.setFont(new Font("Arial", Font.BOLD, 13));
		btnWgtsFile.setEnabled(false);
		btnWgtsFile.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				ChooseFile(txtWgtsfile);
			}
		});
		panel.add(btnWgtsFile, "cell 0 3,growx");
		
		txtWgtsfile = new JTextField();
		txtWgtsfile.setFont(new Font("Arial", Font.PLAIN, 13));
		txtWgtsfile.setEnabled(false);
		txtWgtsfile.setText("weightfile");
		panel.add(txtWgtsfile, "cell 1 3 3 1,growx");
		
		btnOutFile = new JButton("Output File");
		btnOutFile.setFont(new Font("Arial", Font.BOLD, 13));
		btnOutFile.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				ChooseFile(txtOutfile);
			}
		});
		panel.add(btnOutFile, "cell 0 4,growx");

		txtOutfile = new JTextField();
		txtOutfile.setFont(new Font("Arial", Font.PLAIN, 13));
		txtOutfile.setText("outfile");
		panel.add(txtOutfile, "cell 1 4 3 1,growx");

		btnOutTree = new JButton("Output Tree");
		btnOutTree.setFont(new Font("Arial", Font.BOLD, 13));
		btnOutTree.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				ChooseFile(txtOuttree);
			}
		});
		panel.add(btnOutTree, "cell 0 5,growx");

		txtOuttree = new JTextField();
		txtOuttree.setFont(new Font("Arial", Font.PLAIN, 13));
		txtOuttree.setText("outtree");
		panel.add(txtOuttree, "cell 1 5 3 1,growx");
				
		lblUseAnc = new JLabel("Use ancestral states:");
		lblUseAnc.setFont(new Font("Arial", Font.BOLD, 13));
		lblUseAnc.setHorizontalAlignment(SwingConstants.TRAILING);
		panel.add(lblUseAnc, "flowx,cell 0 6 2 1,alignx right");
		
		rdbtnAncY = new JRadioButton("Yes");
		rdbtnAncY.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnAncY.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				AncToggle(true);				
			}
		});
		panel.add(rdbtnAncY, "cell 2 6");
		
		rdbtnAncN = new JRadioButton("No");
		rdbtnAncN.setSelected(true);
		rdbtnAncN.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnAncN.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				AncToggle(false);				
			}
		});
		panel.add(rdbtnAncN, "cell 2 6");
		
		lblUseFac = new JLabel("Use factors:");
		lblUseFac.setFont(new Font("Arial", Font.BOLD, 13));
		lblUseFac.setHorizontalAlignment(SwingConstants.TRAILING);
		panel.add(lblUseFac, "flowx,cell 0 7 2 1,alignx right");
		
		rdbtnFacY = new JRadioButton("Yes");
		rdbtnFacY.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnFacY.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				FacToggle(true);				
			}
		});
		panel.add(rdbtnFacY, "cell 2 7");
		
		rdbtnFacN = new JRadioButton("No");
		rdbtnFacN.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnFacN.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				FacToggle(false);				
			}
		});
		rdbtnFacN.setSelected(true);
		panel.add(rdbtnFacN, "cell 2 7");
		
		lblSiteWgt = new JLabel("Sites weighted:");
		lblSiteWgt.setFont(new Font("Arial", Font.BOLD, 13));
		lblSiteWgt.setHorizontalAlignment(SwingConstants.TRAILING);
		panel.add(lblSiteWgt, "flowx,cell 0 8 2 1,alignx right");
		
		rdbtnWgtY = new JRadioButton("Yes");
		rdbtnWgtY.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnWgtY.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				WgtToggle(true);				
			}
		});
		panel.add(rdbtnWgtY, "cell 2 8");
		
		rdbtnWgtN = new JRadioButton("No");
		rdbtnWgtN.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnWgtN.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				WgtToggle(false);				
			}
		});
		rdbtnWgtN.setSelected(true);
		panel.add(rdbtnWgtN, "cell 2 8");
		
		lblMinCliq = new JLabel("Minimum clique size:");
		lblMinCliq.setFont(new Font("Arial", Font.BOLD, 13));
		lblMinCliq.setHorizontalAlignment(SwingConstants.TRAILING);
		panel.add(lblMinCliq, "flowx,cell 0 9 2 1,alignx right");
		
		rdbtnCliq = new JRadioButton("No");
		rdbtnCliq.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnCliq.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				CliqSizeToggle(rdbtnCliq.isSelected());
			}
		});
		rdbtnCliq.setSelected(true);
		panel.add(rdbtnCliq, "cell 2 9,alignx right");
		
		lblSize = new JLabel("Size:");
		lblSize.setFont(new Font("Arial", Font.BOLD, 13));
		lblSize.setEnabled(false);
		panel.add(lblSize, "cell 3 9");
		
		txtCliqueSize = new JTextField();
		txtCliqueSize.setColumns(5);
		txtCliqueSize.setFont(new Font("Arial", Font.PLAIN, 13));
		txtCliqueSize.setEnabled(false);
		txtCliqueSize.setText("1");
		panel.add(txtCliqueSize, "cell 3 9");
		
		lblOutRoot = new JLabel("Outgroup root:");
		lblOutRoot.setFont(new Font("Arial", Font.BOLD, 13));
		lblOutRoot.setHorizontalAlignment(SwingConstants.TRAILING);
		panel.add(lblOutRoot, "flowx,cell 0 10 2 1,alignx right");
		
		rdbtnOgrpY = new JRadioButton("Yes");
		rdbtnOgrpY.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnOgrpY.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				OutGroupToggle(true);
			}
		});
		panel.add(rdbtnOgrpY, "cell 2 10");
		
		rdbtnOgrpN = new JRadioButton("No");
		rdbtnOgrpN.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnOgrpN.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				OutGroupToggle(false);
			}
		});
		rdbtnOgrpN.setSelected(true);
		panel.add(rdbtnOgrpN, "cell 2 10");
		
		lblOutGroup = new JLabel("Outgroup species:");
		lblOutGroup.setFont(new Font("Arial", Font.BOLD, 13));
		lblOutGroup.setHorizontalAlignment(SwingConstants.TRAILING);
		panel.add(lblOutGroup, "cell 3 10,alignx right");
		
		txtOutRoot = new JTextField();
		txtOutRoot.setColumns(5);
		txtOutRoot.setFont(new Font("Arial", Font.PLAIN, 13));
		txtOutRoot.setText("1");
		panel.add(txtOutRoot, "cell 3 10");
		
		lblMultipleDatasets = new JLabel("Multiple sets in Input file:");
		lblMultipleDatasets.setFont(new Font("Arial", Font.BOLD, 13));
		lblMultipleDatasets.setHorizontalAlignment(SwingConstants.TRAILING);
		panel.add(lblMultipleDatasets, "flowx,cell 0 11 2 1,alignx right");
		
		rdbtnMdsN = new JRadioButton("No");
		rdbtnMdsN.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnMdsN.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				MultiSetToggle(Wgttype.NONE);
			}
		});
		rdbtnMdsN.setSelected(true);
		panel.add(rdbtnMdsN, "cell 2 11,alignx right");
		
		rdbtnMdsY = new JRadioButton("Data sets");
		rdbtnMdsY.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnMdsY.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				MultiSetToggle(Wgttype.DATA);
			}
		});
		panel.add(rdbtnMdsY, "cell 3 11");
		
		rdbtnMwsY = new JRadioButton("Weight sets");
		rdbtnMwsY.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnMwsY.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				MultiSetToggle(Wgttype.WEIGHT);
			}
		});
		panel.add(rdbtnMwsY, "cell 3 11");

		lblNumSets = new JLabel("Number of sets:");
		lblNumSets.setFont(new Font("Arial", Font.BOLD, 13));
		lblNumSets.setHorizontalAlignment(SwingConstants.RIGHT);
		lblNumSets.setEnabled(false);
		panel.add(lblNumSets, "flowx,cell 3 12");
		
		txtSets = new JTextField();
		txtSets.setColumns(5);
		txtSets.setFont(new Font("Arial", Font.PLAIN, 13));
		txtSets.setText("1");
		txtSets.setEnabled(false);
		panel.add(txtSets, "cell 3 12");
		
		lblPrintBefore = new JLabel("Print data before run:");
		lblPrintBefore.setFont(new Font("Arial", Font.BOLD, 13));
		lblPrintBefore.setHorizontalAlignment(SwingConstants.TRAILING);
		panel.add(lblPrintBefore, "flowx,cell 0 13 2 1,alignx right");
		
		rdbtnBrunY = new JRadioButton("Yes");
		rdbtnBrunY.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnBrunY.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				BrunToggle(true);
			}
		});
		panel.add(rdbtnBrunY, "cell 2 13");
		
		rdbtnBrunN = new JRadioButton("No");
		rdbtnBrunN.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnBrunN.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				BrunToggle(false);
			}
		});
		rdbtnBrunN.setSelected(true);
		panel.add(rdbtnBrunN, "cell 2 13");
		
		lblPrintMatrix = new JLabel("Print compatibility matrix:");
		lblPrintMatrix.setFont(new Font("Arial", Font.BOLD, 13));
		lblPrintMatrix.setHorizontalAlignment(SwingConstants.TRAILING);
		panel.add(lblPrintMatrix, "flowx,cell 0 14 2 1,alignx right");
		
		rdbtnCmatY = new JRadioButton("Yes");
		rdbtnCmatY.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnCmatY.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				CmatToggle(true);
			}
		});
		panel.add(rdbtnCmatY, "cell 2 14");
		
		rdbtnCmatN = new JRadioButton("No");
		rdbtnCmatN.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnCmatN.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				CmatToggle(false);
			}
		});
		rdbtnCmatN.setSelected(true);
		panel.add(rdbtnCmatN, "cell 2 14");
		
		lblPrintTree = new JLabel("Print tree:");
		lblPrintTree.setFont(new Font("Arial", Font.BOLD, 13));
		lblPrintTree.setHorizontalAlignment(SwingConstants.TRAILING);
		panel.add(lblPrintTree, "flowx,cell 0 15 2 1,alignx right");
		
		rdbtnPtreeY = new JRadioButton("Yes");
		rdbtnPtreeY.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnPtreeY.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				PtreeToggle(true);
			}
		});
		rdbtnPtreeY.setSelected(true);
		panel.add(rdbtnPtreeY, "cell 2 15");
		
		rdbtnPtreeN = new JRadioButton("No");
		rdbtnPtreeN.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnPtreeN.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				PtreeToggle(false);
			}
		});
		panel.add(rdbtnPtreeN, "cell 2 15");
		
		lblWriteTree = new JLabel("Write tree to Output tree:");
		lblWriteTree.setFont(new Font("Arial", Font.BOLD, 13));
		lblWriteTree.setHorizontalAlignment(SwingConstants.TRAILING);
		panel.add(lblWriteTree, "flowx,cell 0 16 2 1,alignx right");
		
		rdbtnOtreeY = new JRadioButton("Yes");
		rdbtnOtreeY.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnOtreeY.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				OtreeToggle(true);
			}
		});
		rdbtnOtreeY.setSelected(true);
		panel.add(rdbtnOtreeY, "cell 2 16");
		
		rdbtnOtreeN = new JRadioButton("No");
		rdbtnOtreeN.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnOtreeN.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				OtreeToggle(false);
			}
		});
		panel.add(rdbtnOtreeN, "cell 2 16");
		
		lblPrintInd = new JLabel("Display progress:");
		lblPrintInd.setHorizontalAlignment(SwingConstants.RIGHT);
		lblPrintInd.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblPrintInd, "flowx,cell 0 17 2 1,alignx right");
		
		rdbtnPrintIndYes = new JRadioButton("Yes");
		rdbtnPrintIndYes.setBackground(new Color(204, 255, 255));
		rdbtnPrintIndYes.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnPrintIndYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				PrintIndToggle(true);
			}
		});
		rdbtnPrintIndYes.setSelected(true);
		panel.add(rdbtnPrintIndYes, "cell 2 17");
		
		rdbtnPrintIndNo = new JRadioButton("No");
		rdbtnPrintIndNo.setBackground(new Color(204, 255, 255));
		rdbtnPrintIndNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnPrintIndNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				PrintIndToggle(false);
			}
		});
		rdbtnPrintIndNo.setSelected(false);
		panel.add(rdbtnPrintIndNo, "cell 2 17");

		
		btnExecute = new JButton("Execute");
		btnExecute.setFont(new Font("Arial", Font.BOLD, 13));
		btnExecute.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				inputvals = new cliqueData();	
				
				inputvals.infile = txtInputFile.getText();
				inputvals.outfile = txtOutfile.getText();
				inputvals.outfileopt = "w";
				inputvals.outtree = txtOuttree.getText();
				inputvals.outtreeopt = "w";
				inputvals.factfile = txtFactfile.getText();
				inputvals.ancfile = txtAncfile.getText();
				inputvals.weightfile = txtWgtsfile.getText();
				inputvals.useAncStates = rdbtnAncY.isSelected();			
				inputvals.useFactors = rdbtnFacY.isSelected();			
				inputvals.useWeights = rdbtnWgtY.isSelected();			
				inputvals.numWeights = Integer.parseInt(txtSets.getText());					
				inputvals.noCliqueSize = rdbtnCliq.isSelected();				
				inputvals.cliqueSize = Integer.parseInt(txtCliqueSize.getText());
				inputvals.useOutGroupRoot = rdbtnOgrpY.isSelected();			
				inputvals.outRoot = Integer.parseInt(txtOutRoot.getText());				
				inputvals.multipleDataSets = rdbtnMdsY.isSelected();
				inputvals.multipleWeightSets = rdbtnMwsY.isSelected();
				inputvals.printRawData = rdbtnBrunY.isSelected();
				inputvals.progressBar = rdbtnPrintIndYes.isSelected();
				inputvals.printMatrix = rdbtnCmatY.isSelected();
				inputvals.printTree = rdbtnPtreeY.isSelected();
				inputvals.outputTree = rdbtnOtreeY.isSelected();
				
				btnExecute.setEnabled(false);	
				String title = "Clique Progress";
				String curdir = System.getProperty("user.dir");
				curdir += "/progress.txt";
				File fl = new File(curdir);
				fl.delete();
				inTitle = title;
				inCurdir = curdir;
	
				if (checkInputVals())
				{
					
			  	    Thread cliqueThread = new Thread() {
						public void run() {
							runCliqueThreads();
						}
			  	    };
			  	    cliqueThread.start();
				}
				btnExecute.setEnabled(true);
			}
		});
		panel.add(btnExecute, "cell 3 18");
		
		btnQuit = new JButton("Quit");
		btnQuit.setFont(new Font("Arial", Font.BOLD, 13));
		btnQuit.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				if(phylipCall)
				{
					frmCliqueControls.dispose();
					return;
				}
				else
				{
					System.exit(0);
				}
			}
		});
		panel.add(btnQuit, "cell 3 18");
		
	}
	
	protected boolean checkInputVals(){
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
		
		if (inputvals.printTree){
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

		if(inputvals.useAncStates){
			if (!test.FileAvailable(inputvals.ancfile, "Ancestral States"))
			{
				return false;
			}
		}
		
		if(inputvals.useFactors){
			if (!test.FileAvailable(inputvals.factfile, "Factors"))
			{
				return false;
			}
		}
		
		if(inputvals.useWeights){
			if (!test.FileAvailable(inputvals.weightfile, "Weights"))
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
     		new TestFileNames().LibraryMissing("Clique");
    		return;
    	}
		try 
		{
	  	    Thread cliqueRunThread = new Thread() {
	  	    	public void run() {
		  	    	  
		  			// at this point we hook into the C code
		  	   		Clique clique = (Clique) Native.loadLibrary("clique", Clique.class);
		  	        clique.clique(
		  	        		inputvals.infile,
		  	        		inputvals.outfile,
		  	        		inputvals.outfileopt,
		  	        		inputvals.outtree,
		  	        		inputvals.outtreeopt,
		  	        		inputvals.factfile,
		  	        		inputvals.ancfile,
		  	          		inputvals.weightfile,
		  	          		inputvals.useAncStates,
		  	        		inputvals.useFactors,
		  	        		inputvals.useWeights,
		  	        		inputvals.numWeights,
		  	        		inputvals.noCliqueSize,
		  	        		inputvals.cliqueSize,
		  	        		inputvals.multipleDataSets,
		  	        		inputvals.multipleWeightSets,
		  	        		inputvals.useOutGroupRoot,
		  	        		inputvals.outRoot,
		  	        		inputvals.printRawData,
		  	        		inputvals.progressBar,
		  	        		inputvals.printMatrix,
		  	        		inputvals.printTree,
		  	        		inputvals.outputTree);		
					}
		  	    };
	  	  	cliqueRunThread.start();
	  	  	if (inputvals.progressBar)
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
	  	  	if (inputvals.progressBar)
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
