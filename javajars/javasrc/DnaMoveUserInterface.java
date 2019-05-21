package phylip;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Component;
import java.awt.Cursor;
import java.awt.Dimension;
import java.awt.EventQueue;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Point;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.geom.Line2D;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.Scanner;

import javax.swing.AbstractAction;
import javax.swing.Action;
import javax.swing.DefaultComboBoxModel;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JPopupMenu;
import javax.swing.JRadioButton;
import javax.swing.JScrollPane;
import javax.swing.JTextField;
import javax.swing.KeyStroke;
import javax.swing.SwingConstants;

import net.miginfocom.swing.MigLayout;

import utilities.DisplayProgress;

public class DnaMoveUserInterface {
	public class DnaMoveData {
		String infile;
		String intree;
		String outfile;
		String outfileopt;
		String weightfile;
		String outtree;
		String outtreeopt;
	}
	private DnaMoveData inputvals;
	private String inTitle;
	private String inCurdir;
	private DisplayProgress dp;
	private String filedir;
	private boolean phylipCall;

	private JFrame frmDnaMoveControls;
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
	private JLabel lblInitialTreeKind;
	private JComboBox cmbxInitialTreeKind;
	private JButton btnBuildTree;
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
	private JLabel lblSitesWeight;
	private JRadioButton rdbtnSitesYes;
	private JRadioButton rdbtnSitesNo;
	private JLabel lblInputSeq;
	private JRadioButton rdbtnInputSeqYes;
	private JRadioButton rdbtnInputSeqNo;
	private JPanel panTreePlot;

	private JScrollPane paneScrollPane;
	private JPanel panCtlPanel;
	
	private JPopupMenu popAction;
	private JMenuItem pntmMoveNode;
	private JMenuItem pntmOutgroup;
	private JMenuItem pntmTranspose;
	private JMenuItem pntmFlip;
	private JMenuItem pntmUndo;

	private JPopupMenu popMove;
	private JMenuItem mntmMoveTo;
	private JMenuItem mntmMoveBefore;
	private JMenuItem mntmNoMove;

	// node parameters
	Integer m_nodenum;
	NewickNode m_parent;
	ArrayList<NewickNode> m_children;
	String m_name;
	Double m_length;
	NodeType m_type;
	Point m_dispPt;
	Point2D.Double m_corPt;	
	Color m_color;
	
	// Newick tree info
	int m_curr_char;
	char[] m_newick;
	NewickNode root;
	NewickNode newNode;
	NewickNode newBase;
	public enum NodeType{ROOT, NODE, TIP, UNDEF}
	int nodeNum;
	boolean lengthFound;
	
	// plot parameters
	Boolean treeRead;
	int numTips; 
	int maxNameLen;
	ArrayList<NewickNode> treeTips;
	int xOffset;
	int rootLen;
    double xScaleFactor;
    double yScaleFactor;
    int nameOffsetX;
    int nameTweakX;
    int nameTweakY;
    ArrayList<BranchSeg> branchSegs;
    BranchSeg newseg;
    double tipYcor;
    Boolean showLens;
    
    // mouse parameters
    Point2D.Float mousePt;
	int bestTip;
	int bestSeg;
	int moveSeg;
	int moveEnd;
	Component mousePick;
	
	// edit control
	Boolean doMove;
	Cursor oldCursor;
    NewickNode editNode;
	NewickNode longestBranch;
	NewickNode nextLongBranch;

	// Undo array
	ArrayList<String> UndoTrees;

	// keyboard shortcuts
	KeyStroke ctlZ;
	Action doUndo;
	
	protected void IntreeToggle(int selected){
		if(selected == 0){
			txtInputTree.setEnabled(false);
			btnInputTree.setEnabled(false);
	   		UndoTrees.clear();
			frmDnaMoveControls.repaint();        					
		}
		else{
			txtInputTree.setEnabled(true);
			btnInputTree.setEnabled(true);
	   		UndoTrees.clear();
			frmDnaMoveControls.repaint();        					
		}
	}
	
	protected void BuildTree(int selected){
		if(selected == 0){
			createTree();
			frmDnaMoveControls.repaint();        					
		}
		else{
			readTree();
			frmDnaMoveControls.repaint();        					
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

	protected void SiteToggle(boolean isSite) {
		if (isSite) {
			rdbtnSitesYes.setSelected(true);
			rdbtnSitesNo.setSelected(false);
			btnWeightFile.setEnabled(true);
			txtWeightFile.setEnabled(true);
		} else {
			rdbtnSitesYes.setSelected(false);
			rdbtnSitesNo.setSelected(true);
			btnWeightFile.setEnabled(false);
			txtWeightFile.setEnabled(false);
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
	
	/**
	 * Launch the application.
	 */
	public static void main(final String[] args) {
		EventQueue.invokeLater(new Runnable() {
			public void run() {
				try {
					DnaMoveUserInterface window = new DnaMoveUserInterface(args);
					window.frmDnaMoveControls.setVisible(true);
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
		});
	}

	protected void ChooseFile(JTextField file) {
		JFileChooser fileChooser = new JFileChooser(filedir);

		int option = fileChooser.showOpenDialog(frmDnaMoveControls.getRootPane());
		if (option == JFileChooser.APPROVE_OPTION) {
			File selectedFile = fileChooser.getSelectedFile();
			filedir = fileChooser.getCurrentDirectory().getAbsolutePath();
			file.setText(selectedFile.getPath());
		}
	}
	
	/**
	 * Create the application.
	 */
	public DnaMoveUserInterface(String[] args) {
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
	
	@SuppressWarnings("serial")
	private void initialize() {
		
		filedir = System.getProperty("user.dir");
		treeRead = false;
		showLens = false;
		doMove = false;
	    branchSegs = new ArrayList<BranchSeg>();
		UndoTrees = new ArrayList<String> ();	
		
		frmDnaMoveControls = new JFrame();
		frmDnaMoveControls.setBackground(new Color(204, 255, 255));
		frmDnaMoveControls.setTitle("Dnamove");
		frmDnaMoveControls.setBounds(100, 50, 1050, 620);
		frmDnaMoveControls.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frmDnaMoveControls.setPreferredSize(new Dimension(frmDnaMoveControls.getBounds().width, frmDnaMoveControls.getBounds().height));
		
		// register shortcuts
	    int shortcut = Toolkit.getDefaultToolkit().getMenuShortcutKeyMask(); // make key definitions platform independent
		
		ctlZ = KeyStroke.getKeyStroke(KeyEvent.VK_Z,shortcut); // Undo
	    
		doUndo    = new UndoLast();

		paneScrollPane = new JScrollPane();
		paneScrollPane.setPreferredSize(frmDnaMoveControls.getPreferredSize());
		frmDnaMoveControls.getContentPane().add(paneScrollPane);
		
		panCtlPanel = new JPanel();
		panCtlPanel.setPreferredSize(frmDnaMoveControls.getPreferredSize());
		paneScrollPane.setViewportView(panCtlPanel);
		panCtlPanel.setLayout(new MigLayout("", "[pref!,grow][30.00,grow][pref!,grow][pref!,grow][pref!,grow]", "[][][][][][][][][][][grow]"));

		btnInputFile = new JButton("Input File");
		btnInputFile.setFont(new Font("Arial", Font.BOLD, 13));
		btnInputFile.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				ChooseFile(txtInputFile);
			}
		});
		panCtlPanel.add(btnInputFile, "cell 0 0,growx");

		txtInputFile = new JTextField();
		txtInputFile.setFont(new Font("Arial", Font.PLAIN, 13));
		txtInputFile.setText("infile");
		panCtlPanel.add(txtInputFile, "cell 1 0 4 1,growx");
		
		btnInputTree = new JButton("Input Tree");
		btnInputTree.setEnabled(false);
		btnInputTree.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				ChooseFile(txtInputTree);
			}
		});	
		btnInputTree.setFont(new Font("Arial", Font.BOLD, 13));
		panCtlPanel.add(btnInputTree, "cell 0 1,growx");
		
		txtInputTree = new JTextField();
		txtInputTree.setEnabled(false);
		txtInputTree.setFont(new Font("Arial", Font.PLAIN, 13));
		txtInputTree.setText("intree");
		txtInputTree.setBounds(159, 36, 870, 20);
		panCtlPanel.add(txtInputTree, "cell 1 1 4 1,growx");
		
		btnWeightFile = new JButton("Weights File");
		btnWeightFile.setEnabled(false);
		btnWeightFile.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				ChooseFile(txtWeightFile);
			}
		});
		btnWeightFile.setFont(new Font("Arial", Font.BOLD, 13));
		panCtlPanel.add(btnWeightFile, "cell 0 2,growx");
		
		txtWeightFile = new JTextField();
		txtWeightFile.setEnabled(false);
		txtWeightFile.setFont(new Font("Arial", Font.PLAIN, 13));
		txtWeightFile.setText("weightfile");
		txtWeightFile.setBounds(159, 61, 870, 20);
		panCtlPanel.add(txtWeightFile, "cell 1 2 4 1,growx");
	
		btnOutputFile = new JButton("Output File");
		btnOutputFile.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				ChooseFile(txtOutputFile);
			}
		});
		btnOutputFile.setFont(new Font("Arial", Font.BOLD, 13));
		panCtlPanel.add(btnOutputFile, "cell 0 3,growx");
		
		txtOutputFile = new JTextField();
		txtOutputFile.setFont(new Font("Arial", Font.PLAIN, 13));
		txtOutputFile.setText("outfile");
		txtOutputFile.setBounds(159, 86, 870, 20);
		panCtlPanel.add(txtOutputFile, "cell 1 3 4 1,growx");
		
		btnOutputTree = new JButton("Output Tree");
		btnOutputTree.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				ChooseFile(txtOutputTree);
			}
		});
		btnOutputTree.setFont(new Font("Arial", Font.BOLD, 13));
		panCtlPanel.add(btnOutputTree, "cell 0 4,growx");
		
		txtOutputTree = new JTextField();
		txtOutputTree.setText("outtree");
		txtOutputTree.setFont(new Font("Arial", Font.PLAIN, 13));
		txtOutputTree.setBounds(159, 111, 870, 20);
		panCtlPanel.add(txtOutputTree, "cell 1 4 4 1,growx");
		
		lblInitialTreeKind = new JLabel("Initial tree:");
		lblInitialTreeKind.setHorizontalAlignment(SwingConstants.RIGHT);
		lblInitialTreeKind.setFont(new Font("Arial", Font.BOLD, 13));
		panCtlPanel.add(lblInitialTreeKind, "cell 0 5 1 1,alignx right");
		
		cmbxInitialTreeKind = new JComboBox();
		cmbxInitialTreeKind.setModel(new DefaultComboBoxModel(new String[] {"Arbitrary", "User trees in input file"}));
		cmbxInitialTreeKind.setSelectedIndex(0);
		cmbxInitialTreeKind.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				IntreeToggle(cmbxInitialTreeKind.getSelectedIndex());
			}
		});
		panCtlPanel.add(cmbxInitialTreeKind, "flowx,cell 1 5");
		
		btnBuildTree = new JButton("Build Tree");
		btnBuildTree.setFont(new Font("Arial", Font.BOLD, 13));
		btnBuildTree.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				BuildTree(cmbxInitialTreeKind.getSelectedIndex());
			}
		});
		panCtlPanel.add(btnBuildTree, "cell 1 5");

		lblInputSeq = new JLabel("Input sequences:");
		lblInputSeq.setEnabled(true);
		lblInputSeq.setFont(new Font("Arial", Font.BOLD, 13));
		lblInputSeq.setHorizontalAlignment(SwingConstants.RIGHT);
		panCtlPanel.add(lblInputSeq, "cell 0 6,alignx right");

		rdbtnInputSeqYes = new JRadioButton("Interleaved");
		rdbtnInputSeqYes.setEnabled(true);
		rdbtnInputSeqYes.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnInputSeqYes.setHorizontalAlignment(SwingConstants.LEFT);
		rdbtnInputSeqYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				InputSeqToggle(true);
			}
		});
		rdbtnInputSeqYes.setSelected(true);
		panCtlPanel.add(rdbtnInputSeqYes, "cell 1 6");
		
		lblSitesWeight = new JLabel("Sites weighted:");
		lblSitesWeight.setFont(new Font("Arial", Font.BOLD, 13));
		lblSitesWeight.setHorizontalAlignment(SwingConstants.RIGHT);
		panCtlPanel.add(lblSitesWeight, "flowx,cell 0 7,alignx right");
	
		rdbtnSitesYes = new JRadioButton("Yes");
		rdbtnSitesYes.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnSitesYes.setHorizontalAlignment(SwingConstants.LEFT);
		rdbtnSitesYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				SiteToggle(true);
			}
		});
		rdbtnSitesYes.setSelected(false);
		panCtlPanel.add(rdbtnSitesYes, "cell 1 7");
	
		rdbtnSitesNo = new JRadioButton("No");
		rdbtnSitesNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnSitesNo.setHorizontalAlignment(SwingConstants.LEFT);
		rdbtnSitesNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				SiteToggle(false);
			}
		});
		rdbtnSitesNo.setSelected(true);
		panCtlPanel.add(rdbtnSitesNo, "cell 1 7");

		lblOutRoot = new JLabel("Outgroup Root:");
		lblOutRoot.setFont(new Font("Arial", Font.BOLD, 13));
		lblOutRoot.setHorizontalAlignment(SwingConstants.RIGHT);
		panCtlPanel.add(lblOutRoot, "flowx,cell 2 5,alignx right");

		rdbtnOutYes = new JRadioButton("Yes");
		rdbtnOutYes.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnOutYes.setHorizontalAlignment(SwingConstants.LEFT);
		rdbtnOutYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				OutToggle(true);
			}
		});
		rdbtnOutYes.setSelected(false);
		panCtlPanel.add(rdbtnOutYes, "cell 3 5");

		rdbtnOutNo = new JRadioButton("No, use as outgroup species");
		rdbtnOutNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnOutNo.setHorizontalAlignment(SwingConstants.LEFT);
		rdbtnOutNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				OutToggle(false);
			}
		});
		rdbtnOutNo.setSelected(true);
		panCtlPanel.add(rdbtnOutNo, "cell 3 5");

		lblOutNumber = new JLabel("Number of the outgroup:\r\n");
		lblOutNumber.setFont(new Font("Arial", Font.BOLD, 13));
		lblOutNumber.setHorizontalAlignment(SwingConstants.RIGHT);
		lblOutNumber.setEnabled(false);
		panCtlPanel.add(lblOutNumber, "flowx,cell 2 6,alignx right");

		txtOutNumber = new JTextField();
		txtOutNumber.setText("1");
		txtOutNumber.setFont(new Font("Arial", Font.PLAIN, 13));
		txtOutNumber.setEnabled(false);
		txtOutNumber.setColumns(6);
		panCtlPanel.add(txtOutNumber, "cell 3 6");

		lblUseThreshold = new JLabel("Use Threshold parsimony:");
		lblUseThreshold.setFont(new Font("Arial", Font.BOLD, 13));
		lblUseThreshold.setHorizontalAlignment(SwingConstants.RIGHT);
		panCtlPanel.add(lblUseThreshold, "flowx,cell 2 7,alignx right");

		rdbtnThresholdYes = new JRadioButton("Yes");
		rdbtnThresholdYes.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnThresholdYes.setHorizontalAlignment(SwingConstants.LEFT);
		rdbtnThresholdYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				ThresholdToggle(true);
			}
		});
		rdbtnThresholdYes.setSelected(false);
		panCtlPanel.add(rdbtnThresholdYes, "cell 3 7");

		rdbtnThresholdNo = new JRadioButton("No, use ordinary parsimony");
		rdbtnThresholdNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnThresholdNo.setHorizontalAlignment(SwingConstants.LEFT);
		rdbtnThresholdNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				ThresholdToggle(false);
			}
		});
		rdbtnThresholdNo.setSelected(true);
		panCtlPanel.add(rdbtnThresholdNo, "cell 3 7");

		lblThresholdValue = new JLabel("Threshold value:");
		lblThresholdValue.setFont(new Font("Arial", Font.BOLD, 13));
		lblThresholdValue.setHorizontalAlignment(SwingConstants.RIGHT);
		lblThresholdValue.setEnabled(false);
		panCtlPanel.add(lblThresholdValue, "flowx,cell 2 8,alignx right");

		txtThresholdValue = new JTextField();
		txtThresholdValue.setText("1.0");
		txtThresholdValue.setFont(new Font("Arial", Font.PLAIN, 13));
		txtThresholdValue.setEnabled(false);
		txtThresholdValue.setColumns(6);
		panCtlPanel.add(txtThresholdValue, "cell 3 8");

		rdbtnInputSeqNo = new JRadioButton("Sequential");
		rdbtnInputSeqNo.setEnabled(true);
		rdbtnInputSeqNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnInputSeqNo.setHorizontalAlignment(SwingConstants.LEFT);
		rdbtnInputSeqNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				InputSeqToggle(false);
			}
		});
		rdbtnInputSeqNo.setSelected(false);
		panCtlPanel.add(rdbtnInputSeqNo, "cell 1 6");
		
		btnExecute = new JButton("Execute");
		btnExecute.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				inputvals = new DnaMoveData();
				
				btnExecute.setEnabled(false);	
				String title = "Dnamove Progress";
				String curdir = System.getProperty("user.dir");
				curdir += "/progress.txt";
				File fl = new File(curdir);
				fl.delete();
				inTitle = title;
				inCurdir = curdir;
	/*
				if (checkInputVals())
				{
					
			  	    Thread dnaParsThread = new Thread() {
						public void run() {
							runDnaParsThreads();
						}
			  	    };
			  	    dnaParsThread.start();
				}
				*/
				btnExecute.setEnabled(true);
			}
		});
		btnExecute.setFont(new Font("Arial", Font.BOLD, 13));
		panCtlPanel.add(btnExecute, "cell 2 9,alignx right");

		btnQuit = new JButton("Quit");
		btnQuit.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				if(phylipCall)
				{
					frmDnaMoveControls.dispose();
					return;
				}
				else
				{
					System.exit(0);
				}
			}
		});
		btnQuit.setFont(new Font("Arial", Font.BOLD, 13));
		panCtlPanel.add(btnQuit, "cell 3 9");
	
 	 	panTreePlot = new JPanel()
        {
   	 		
            public void paintComponent(Graphics graph)
            {
                draw(graph);
            }
            
        };
	    panTreePlot.setBackground(Color.WHITE); // does not seem to work JRMFix
		panCtlPanel.add(panTreePlot, "cell 0 10 5 1,grow");
		
		oldCursor = panTreePlot.getCursor();
		panTreePlot.addMouseListener(new MouseAdapter() {
			public void mouseClicked(MouseEvent e){
				if(doMove)
				{
					popMove.show(e.getComponent(), e.getX(), e.getY());
					mousePt = new Point2D.Float(e.getX(), e.getY());
					mousePick = e.getComponent();
					panTreePlot.setCursor(oldCursor);
					double minsegdist = 100000;
					moveSeg = 0;
					for(int i=0; i<branchSegs.size(); i++)
					{
						BranchSeg seg = branchSegs.get(i);
						if ((seg.distanceToBranch(mousePt)) < minsegdist)
						{
							minsegdist = seg.distanceToBranch(mousePt);
							moveSeg = i;
						}
					}		
			        branchSegs.get(moveSeg).m_endnode.m_color = Color.GREEN;
					doMove = false;
					frmDnaMoveControls.repaint();        					
				}
				else
				{
					if (root != null)
						{
		 				root.ClearTree(root);
						popAction.show(e.getComponent(), e.getX(), e.getY());
						mousePt = new Point2D.Float(e.getX(), e.getY());
						//System.out.println("Mouse Clicked at X: " + e.getX() + " Y: " + e.getY());
						double mintipdist = 100000; 
						bestTip = 0;
						for(int i=0; i<treeTips.size(); i++)
						{
							NewickNode node = treeTips.get(i);
							if (node.m_dispPt.distance(mousePt) < mintipdist)
							{
								mintipdist = node.m_dispPt.distance(mousePt);
								bestTip = i;
							}
						}
						
						double minsegdist = 100000;
						bestSeg = 0;
						for(int i=0; i<branchSegs.size(); i++)
						{
							BranchSeg seg = branchSegs.get(i);
							if ((seg.distanceToBranch(mousePt)) < minsegdist)
							{
								minsegdist = seg.distanceToBranch(mousePt);
								bestSeg = i;
							}
						}
						
						if (mintipdist <= minsegdist)
						{
							pntmFlip.setEnabled(false);
							pntmTranspose.setEnabled(false);
					        pntmMoveNode.setEnabled(false);
					        pntmOutgroup.setEnabled(false);
					        treeTips.get(bestTip).m_color = Color.BLUE;
						}
						else
						{
					        pntmOutgroup.setEnabled(true);
					        
					        branchSegs.get(bestSeg).m_endnode.m_color = Color.BLUE;
					        if (branchSegs.get(bestSeg).m_endnode.m_type != NodeType.TIP)
					        {
								pntmFlip.setEnabled(true);
								pntmTranspose.setEnabled(true);
					        }
					        else
					        {
								pntmFlip.setEnabled(false);
								pntmTranspose.setEnabled(false);
					        }
					        
					        if (branchSegs.get(bestSeg).m_endnode.m_type == NodeType.ROOT)
					        {
						        pntmMoveNode.setEnabled(false);				        	
					        }
					        else
					        {
						        pntmMoveNode.setEnabled(true);
					        }
						}
						frmDnaMoveControls.repaint();        					
					}	
				}
			}
	    });
	    
        popAction= new JPopupMenu();
        
        pntmMoveNode = new JMenuItem("Move node");
        pntmMoveNode.setEnabled(false);
        pntmMoveNode.addMouseListener(new MouseAdapter() {
			public void mouseClicked(MouseEvent e){
 				root.ClearTree(root);
				frmDnaMoveControls.repaint();
			}
		});
        pntmMoveNode.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				doMove = true;
				panTreePlot.setCursor(Cursor.getPredefinedCursor(Cursor.CROSSHAIR_CURSOR));
		    }
		});
        
        popAction.add(pntmMoveNode);

        popMove= new JPopupMenu();

        mntmMoveBefore = new JMenuItem("Before node");
        mntmMoveBefore.setEnabled(true);
        mntmMoveBefore.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
	       		NewickNode movingNode = branchSegs.get(bestSeg).m_endnode;
	       		NewickNode moveBeforeNode = branchSegs.get(moveSeg).m_endnode;
	       		moveNodeBefore(movingNode, moveBeforeNode);
				Redisplay(true);
		    }
		});        
        popMove.add(mntmMoveBefore);

        mntmMoveTo = new JMenuItem("To node");
        mntmMoveTo.setEnabled(true);
        mntmMoveTo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
	       		NewickNode movingNode = branchSegs.get(bestSeg).m_endnode;
	       		NewickNode moveToNode;
	       		if ((Math.abs(branchSegs.get(moveSeg).m_startnode.m_dispPt.x-mousePt.x)) < 
	       			 Math.abs((branchSegs.get(moveSeg).m_endnode.m_dispPt.x-mousePt.x)))
       			{
	       			moveToNode = branchSegs.get(moveSeg).m_startnode;
       			}
	       		else if(branchSegs.get(moveSeg).m_endnode.m_type == NodeType.TIP)
	       		{
	       			moveToNode = branchSegs.get(moveSeg).m_startnode;	       			
	       		}
	       		else
	       		{
	       			moveToNode = branchSegs.get(moveSeg).m_endnode;	       				       			
	       		}
	       		moveNodeTo(movingNode, moveToNode);
  				calcNodeCoords(root);
				Redisplay(true);
		    }
		});        
        popMove.add(mntmMoveTo);

        mntmNoMove = new JMenuItem("Cancel");
        mntmNoMove.setEnabled(true);
        mntmNoMove.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				root.ClearTree(root);
				frmDnaMoveControls.repaint();
		    }
		});        
        popMove.add(mntmNoMove);

        pntmOutgroup = new JMenuItem("Outgroup");
        pntmOutgroup.setEnabled(false);
        pntmOutgroup.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				doMove = false;
	       		NewickNode movingNode = branchSegs.get(bestSeg).m_endnode;
	       		outgroupRoot(movingNode);
				Redisplay(true);
		    }
		});
        popAction.add(pntmOutgroup);

        pntmTranspose = new JMenuItem("Transpose");
        pntmTranspose.setEnabled(false);
        pntmTranspose.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
	       		editNode = branchSegs.get(bestSeg).m_endnode;
	       		transposeBranch(editNode);
				Redisplay(true);
		    }
		});
        popAction.add(pntmTranspose);

        pntmFlip = new JMenuItem("Flip");
        pntmFlip.setEnabled(false);
        pntmFlip.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
	       		editNode = branchSegs.get(bestSeg).m_endnode;
	       		flipBranch(editNode);
				Redisplay(true);
		    }
		});
        
        popAction.add(pntmFlip);

        pntmUndo = new JMenuItem("Undo");
        pntmUndo.setEnabled(false);
        pntmUndo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				undoLast();
		    }
		});
        popAction.add(pntmUndo);     

	}
	
	public void readTree() {
		String intree = txtInputTree.getText();
		String newickTree = "";
		nodeNum = -1;
	    try 
	    {
			Scanner sf = new Scanner(new File(intree));
			while (sf.hasNextLine()) 
			{
				String curline = sf.nextLine();
				newickTree += curline;
			}
			parseNewick(newickTree);
	    }
		catch (FileNotFoundException e)
		{
			String msg = "Input tree: ";
			msg += intree;
			msg += " does not exist.";
			JOptionPane.showMessageDialog(null, msg, "Error", JOptionPane.ERROR_MESSAGE);			
		}
   		UndoTrees.clear();
		addToUndo(newickTree);
	}
	
	public void createTree() {
		String intree = txtInputFile.getText();
		nodeNum = -1;
		ArrayList<String> names = new ArrayList<String> ();	
		String newickTree = "";
	    try 
	    {
			Scanner sf = new Scanner(new File(intree));
			int nline = 0;
			while (sf.hasNextLine()) 
			{
	        	Scanner linescan =  new Scanner( sf.nextLine());
	        	String inline = linescan.next();
				if (nline > 0)
				{
		        	names.add(inline);
				}
				nline++;
			}
			newickTree = buildNewick(names);
			parseNewick(newickTree);
	    }
		catch (FileNotFoundException e)
		{
			String msg = "Input tree: ";
			msg += intree;
			msg += " does not exist.";
			JOptionPane.showMessageDialog(null, msg, "Error", JOptionPane.ERROR_MESSAGE);			
		}
   		UndoTrees.clear();
		addToUndo(newickTree);
	}
	
	public String buildNewick(ArrayList<String> names)
	{
		String newTree = "";
		for (int i=0; i<names.size(); i++)
		{
			if (i != names.size()-1)
			{
				newTree +="(";
				newTree += names.get(i);
				newTree += ",";
			}
			else
			{
				newTree += names.get(i);
			}		
		}
		for (int i=0; i<names.size()-1; i++)
		{
			newTree += ")";
		}		
		newTree += ";";
		return newTree;
	}
	public void calcNodeCoords(NewickNode root)
	{
		// first fill in xcor
		walkDownFromRoot(root);
		
		// now do ycor
		for(int i=0; i<treeTips.size(); i++)
		{
			walkUpFromTip(treeTips.get(i));
		}
	}
	
	public void walkDownFromRoot(NewickNode node)
	{
		// calculates x coordinates
		node.m_corPt.x = 0.0;
		NewickNode parent = node.GetParent();
		if (parent != null)
		{
			node.m_corPt.x = parent.m_corPt.x + node.m_length;
		}

		if (node.m_type == NodeType.TIP)
		{
			return;
		}
		else
		{
			Iterator<NewickNode> itr = node.m_children.iterator();
			while(itr.hasNext())
			{
				walkDownFromRoot(itr.next());
			}
		}
	}
	
	public void walkUpFromTip(NewickNode tip)
	{
		// calculates y coordinates
		NewickNode parent = tip.GetParent();
		if (parent != null)
		{
			if (parent.m_corPt.y != 0.0)
			{
				walkUpFromTip(parent);
			}
			else
			{
				walkDownFromNode(parent);
			}
		}
	}
	
	public void walkDownFromNode(NewickNode node)
	{
		Iterator<NewickNode> itr = node.m_children.iterator();
		double ySum = 0.0;
		while(itr.hasNext())
		{
			double ynew = findLowestChildY(itr.next());
			ySum += ynew;
		}
		node.m_corPt.y = ySum/node.m_children.size();
	}
	
	public double findLowestChildY(NewickNode node)
	{
		if (node.m_type != NodeType.TIP)
		{
			if (node.m_corPt.y == 0.0)
			{
				walkDownFromNode(node);
			}
		}
		return node.m_corPt.y;
	}
	
	public void moveNodeBefore(NewickNode movingNode, NewickNode moveBeforeNode)
	{
		if ((movingNode.m_parent == moveBeforeNode.m_parent) && (movingNode.m_parent.m_children.size() == 2))
		{
			String msg = "Moving a node before a node that shares the moving node's parent is not meaningful.";
			JOptionPane.showMessageDialog(null, msg, "Error", JOptionPane.ERROR_MESSAGE);

			return;
		}
		if ((movingNode.m_parent == root) && (root.m_children.size() == 2))
		{	
			// root case with only 2 children is different because the root moves
			NewickNode newRoot = new NewickNode();
			Iterator<NewickNode> itr = root.m_children.iterator();
			while(itr.hasNext())
			{
				NewickNode child = itr.next();
				if (child != movingNode)
				{
					movingNode.m_length += child.m_length;
					newRoot = child;
					
					// clear things out
					newRoot.m_length = 1.0;
					newRoot.m_type = NodeType.ROOT;
					newRoot.m_parent = null;
					newRoot.m_dispPt = new Point(0,0);
					newRoot.m_corPt = new Point2D.Double(0.0, 0.0);
				}
			}

			NewickNode newNode = moveBeforeNode.m_parent.AddChild();
			movingNode.m_parent = newNode;
			disconnectNode(moveBeforeNode);
			moveBeforeNode.m_parent = newNode;
			moveBeforeNode.m_length /= 2;
			newNode.m_length = moveBeforeNode.m_length;
			newNode.m_corPt.x = newNode.m_parent.m_corPt.x + newNode.m_length;
			newNode.m_children.add(moveBeforeNode);
			newNode.m_children.add(movingNode);
			root  = newRoot;
		}
		else
		{
			NewickNode newNode = moveBeforeNode.m_parent.AddChild();
			disconnectNode(movingNode);
			movingNode.m_parent = newNode;
			disconnectNode(moveBeforeNode);
			moveBeforeNode.m_parent = newNode;
			moveBeforeNode.m_length /= 2;
			newNode.m_length = moveBeforeNode.m_length;
			newNode.m_corPt.x = newNode.m_parent.m_corPt.x + newNode.m_length;
			newNode.m_children.add(moveBeforeNode);
			newNode.m_children.add(movingNode);
		}
		rebuildTipList();
		recalcYcors();
		calcNodeCoords(root);
	}
	
	public void moveNodeTo(NewickNode movingNode, NewickNode moveToNode)
	{
		if (movingNode.m_parent == moveToNode)
		{
			String msg = "Moving a node to it's own parent is not meaningful.";
			JOptionPane.showMessageDialog(null, msg, "Error", JOptionPane.ERROR_MESSAGE);

			return;
		}
		
		if (movingNode.m_parent == root) 
		{	
			// root case with only 2 children is different because the root moves
			NewickNode newRoot = new NewickNode();
			Iterator<NewickNode> itr = root.m_children.iterator();
			while(itr.hasNext())
			{
				NewickNode child = itr.next();
				if (child != movingNode)
				{
					movingNode.m_length += child.m_length;
					newRoot = child;
					
					// clear things out
					newRoot.m_length = 1.0;
					newRoot.m_type = NodeType.ROOT;
					newRoot.m_parent = null;
					newRoot.m_dispPt = new Point(0,0);
					newRoot.m_corPt = new Point2D.Double(0.0, 0.0);
				}
			}
			disconnectNode(movingNode);
			movingNode.m_parent = moveToNode;
			moveToNode.m_children.add(movingNode);
			root  = newRoot;
		}
		else
		{
			disconnectNode(movingNode);
			movingNode.m_parent = moveToNode;
			moveToNode.m_children.add(movingNode);
		}
		rebuildTipList();
		recalcYcors();
		calcNodeCoords(root);
	}
	
	public void disconnectNode(NewickNode movingNode)
	{
		if (movingNode.m_parent.m_children.size() > 2)
		{
			// remove child
			Iterator<NewickNode> itr = movingNode.m_parent.m_children.iterator();
			while(itr.hasNext())
			{
				if (itr.next() == movingNode)
				{
					itr.remove();
				}
			}
		}
		else if (movingNode.m_parent.m_children.size() == 2)
		{
			if (movingNode.m_parent.m_type == NodeType.ROOT)
			{
				Iterator<NewickNode> itr = movingNode.m_parent.m_children.iterator();
				while(itr.hasNext())
				{
					NewickNode child = itr.next();
					if (child == movingNode)
					{
						itr.remove();
					}
				}		
			}
			else
			{
				// removing one child also eliminates the parent
				// transfer the length and the not moving child to the grandparent
				Iterator<NewickNode> itr = movingNode.m_parent.m_children.iterator();
				while(itr.hasNext())
				{
					NewickNode child = itr.next();
					if (child != movingNode)
					{
						// if .m_parent.m_parent == null, at the root and the rules change
						child.m_length += movingNode.m_parent.m_length;
						child.m_parent = movingNode.m_parent.m_parent;
						movingNode.m_parent.m_parent.m_children.add(child);
					}
				}	
				
				// eliminate the parent
				itr = movingNode.m_parent.m_parent.m_children.iterator();
				while(itr.hasNext())
				{
					NewickNode child = itr.next();
					if (child == movingNode.m_parent)
					{
						itr.remove();
						break;
					}
				}
			}
		}
		else
		{
			System.out.println("Somehow parent " + movingNode.m_parent.m_nodenum + " slipped through with a single child.");
		}		
	}

	public void transposeBranch(NewickNode node)
	{
		Collections.reverse(node.m_children);
		recalcYcors();
	}
	
	public void recalcYcors()
	{
		tipYcor = .5;
		resetYcors(root);
		for(int i=0; i<treeTips.size(); i++)
		{
			walkUpFromTip(treeTips.get(i));
		}		
	}
	
	public void resetYcors(NewickNode node)
	{
		// this prepares the tree so the internal nodes will be recalculated
		if (node.m_type == NodeType.TIP)
		{
			node.m_corPt.y = tipYcor;
			tipYcor += 1.0;
			return;
		}
		else
		{
			node.m_corPt.y = 0.0;
			Iterator<NewickNode> itr = node.m_children.iterator();
			while(itr.hasNext())
			{
				resetYcors(itr.next());
			}
		}		
	}
		
	void draw(Graphics g) {
	    Graphics2D g2d = (Graphics2D) g;
		g2d.setColor(frmDnaMoveControls.getBackground());
		g2d.fillRect(0, 0, frmDnaMoveControls.getWidth(), frmDnaMoveControls.getHeight());
		g2d.setColor(Color.BLACK);

		plotDnaTree(g2d);
	}
	void plotDnaTree(Graphics2D g2d){
		//BasicStroke stroke1 = new BasicStroke(1.0f, BasicStroke.CAP_BUTT, BasicStroke.JOIN_MITER, 10.0f, new float[] {4.0f, 4.0f}, 0.0f);
		//BasicStroke stroke2 = new BasicStroke(1.0f, BasicStroke.CAP_BUTT, BasicStroke.JOIN_MITER, 1.0f);
		BasicStroke stroke3 = new BasicStroke(5.0f, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND, 3.0f);
	    Font plotFont = new Font("helvetica", Font.BOLD, 14);
	    g2d.setFont(plotFont);
		g2d.setStroke(stroke3);
		//g2d.drawLine(0, 0, 100, 100);
	    if (treeRead)
	    {
	    	xOffset = 20;
	    	rootLen = 30;
	    	nameOffsetX = 5;
	    	nameTweakX = 30;
	    	nameTweakY = 5;
		    			
			// find longest path and longest label
			FontMetrics fm = g2d.getFontMetrics(plotFont);
			int nameWid = 0;
			double pathWid = 0.0;
			for(int i=0; i<treeTips.size(); i++)
			{
				if (fm.stringWidth(treeTips.get(i).m_name) > nameWid)
				{
					nameWid = fm.stringWidth(treeTips.get(i).m_name);
				}
				if (treeTips.get(i).m_corPt.x > pathWid)
				{
					pathWid = treeTips.get(i).m_corPt.x;
				}
			}
			
		    xScaleFactor = (panTreePlot.getWidth() - (nameWid + xOffset + rootLen + nameTweakX)) / pathWid;
		    yScaleFactor = panTreePlot.getHeight() / (double)(treeTips.size() + 1);
						
			// plot tree
		    if (!branchSegs.isEmpty())
		    {
		    	branchSegs.clear();
		    }
			root.PlotTree(root, g2d);
	    }
	}
	
	public void parseNewick(String tree){
		m_newick = new char[tree.length()];
		tree.getChars(0, tree.length(), m_newick, 0);
		m_curr_char = 0;
		root = new NewickNode();
		root.m_type = NodeType.ROOT;
	    nodeNum = 0;
	    root.m_nodenum = nodeNum;
		NewickNode current = root;
	    treeTips = new ArrayList<NewickNode>();	
	    lengthFound = false; // if even one is found on input, the whole output tree gets lengths
		
		// parameters needed for plot
		numTips = 0;
		maxNameLen = 0;
		treeRead = false;
	
		while (m_newick[m_curr_char] != ';')
		{
			switch(m_newick[m_curr_char])
			{
            case '(':
                ++m_curr_char;
                // create first daughter
                current = current.AddChild();
                //System.out.print("found first daughter: ");
                break;
            case ',':
                ++m_curr_char;
                // create additional daughter
                current = current.GetParent().AddChild();
                //System.out.print("found second daughter: ");
                break;
            case ')':
                ++m_curr_char;
                // coalesce daughters
                current = current.GetParent().Terminate();
                //System.out.println("end of pair");
                break;
            case ' ':
            case '\n':
            case '\t':
            case '\r':
                ++m_curr_char;
                // skip whitespace
                break;
            case ':':
                current.ProcessLength();
                break;
            case '[':
                 current.ProcessComment();
                 break;
            default:
                // Anything unrecognized must be a tip name
                current.ProcessName();
                current.m_type = NodeType.TIP;
                current.m_corPt.y = (double)numTips + .5;
                treeTips.add(current);
                numTips++;
				}	
		}
		calcNodeCoords(root);
		//root.DumpTree(root);
		treeRead = true;
	}
	
	
	public void undoLast()
	{
		if (UndoTrees.size() > 1)
		{
			UndoTrees.remove(UndoTrees.size()-1);
			parseNewick(UndoTrees.get(UndoTrees.size()-1));
			calcNodeCoords(root);
			frmDnaMoveControls.repaint();
			if(UndoTrees.size() == 1)
			{
		        pntmUndo.setEnabled(false);						
			}
		}		
	}
	
	@SuppressWarnings("serial")
	public class UndoLast extends AbstractAction{
		public UndoLast(){
		}
		public void actionPerformed(ActionEvent e){
			undoLast();
			doMove = false;
		}
	};
	
	public void rebuildTipList()
	{
		treeTips.clear();
		Iterator<NewickNode> itr = root.m_children.iterator();
		while(itr.hasNext())
		{
			NewickNode child = itr.next();
			findNextTip(child);
		}		
	}
	
	public void findNextTip(NewickNode node)
	{
		if (node.m_type == NodeType.TIP)
		{
			treeTips.add(node);
			return;
		}
		else
		{
			Iterator<NewickNode> itr = node.m_children.iterator();
			while(itr.hasNext())
			{
				findNextTip(itr.next());
			}		
		}
	}

	public void flipBranch(NewickNode node)
	{
		Collections.reverse(node.m_children);
		Iterator<NewickNode> itr = node.m_children.iterator();
		while(itr.hasNext())
		{
			NewickNode child = itr.next();
			flipBranch(child);
		}		
		recalcYcors();
	}
	
	public void addToUndo(String newTree)
	{
		UndoTrees.add(newTree);
		if (UndoTrees.size() > 1)
		{
			pntmUndo.setEnabled(true);
		}
		else
		{
	        pntmUndo.setEnabled(false);
		}
	}
	
	public void outgroupRoot(NewickNode movingNode)
	{
		disconnectNode(movingNode);
		movingNode.m_length/=2;
		NewickNode newRoot = new NewickNode();
		//newRoot.NewRoot(newRoot);
		nodeNum++;
		newRoot.m_nodenum = nodeNum;
		newRoot.m_type = NodeType.ROOT;
		movingNode.m_parent = newRoot;
		
		NewickNode oldRoot = root;
		oldRoot.m_type = NodeType.NODE;
		oldRoot.m_parent = newRoot;
		oldRoot.m_length = movingNode.m_length;

		newRoot.m_children.add(oldRoot);
		newRoot.m_children.add(movingNode);
		root = newRoot;
		
		recalcYcors();
		calcNodeCoords(root);
	}
	
	public void Redisplay(boolean doUndo)
	{
		if (doUndo)
		{
			addToUndo(toNewick());
		}
		root.ClearTree(root);
		frmDnaMoveControls.repaint();
	}
	
	public String toNewick()
	{
		String retstr = "";
		retstr = addToNewick(root, retstr);
		
		return retstr;
	}
	
	public String addToNewick(NewickNode node, String retstr)
	{
		if(node.m_type == NodeType.TIP)
		{
			retstr += node.m_name;
			if (lengthFound)
			{
				retstr += ":"+node.m_length;
			}
		}
		else if (node.m_type == NodeType.ROOT)
		{
			Iterator<NewickNode> itr = node.m_children.iterator();
			retstr += "(";
			int numnode = 0;
			while(itr.hasNext())
			{
				if (numnode > 0)
				{
					retstr += ",";
				}
				retstr = addToNewick(itr.next(), retstr);
				numnode++;
			}
			retstr += ");\n";
		}
		else //node.m_type == NodeType.NODE
		{
			Iterator<NewickNode> itr = node.m_children.iterator();
			retstr += "(";
			int numnode = 0;
			while(itr.hasNext())
			{
				if (numnode > 0)
				{
					retstr += ",";
				}
				retstr = addToNewick(itr.next(), retstr);
				numnode++;
			}
			retstr += ")";
			if (lengthFound)
			{
				retstr += ":"+node.m_length;
			}
		}
		return retstr;
	}
	
	public class NewickNode{
		Integer m_nodenum;
		NewickNode m_parent;
		ArrayList<NewickNode> m_children;
		String m_name;
		Double m_length;
		NodeType m_type;
		Point m_dispPt;
		Point2D.Double m_corPt;	
		Color m_color;
		
		public NewickNode()
		{
			m_nodenum = -1;
			m_parent = null;
		    m_children = new ArrayList<NewickNode>();			
			m_name = "";
			m_length = 1.0;
			m_type = NodeType.UNDEF;
			m_dispPt = new Point(0,0);
			m_corPt = new Point2D.Double(0.0, 0.0);
			m_color = Color.BLACK;
		}
		
		public NewickNode AddChild() {
		    newNode = new NewickNode();
		    m_children.add(newNode);
		    nodeNum++;
		    newNode.m_nodenum = nodeNum;
		    newNode.m_parent = this;
		    newNode.m_type = NodeType.NODE;
		    return newNode;
		}
		
		public NewickNode NewRoot(NewickNode node) {
			nodeNum++;
			node.m_nodenum = nodeNum;
			node.m_type = NodeType.ROOT;
		    return node;
		}

		public NewickNode GetParent(){
			return m_parent;
		} 
		
		public NewickNode Terminate()
		{
		    assert(m_children.size() > 1);
		    return this;
		} 

		public void ProcessLength()
		{
			++m_curr_char;  // skip the colon
			
			String lenval = "";
			int lenstr = 0;
            while ((m_newick[m_curr_char] != ',') &&
             	   (m_newick[m_curr_char] != ')'))
            {
            	lenval += m_newick[m_curr_char];
            	++m_curr_char;
            	lenstr++;
            }
            if (lenstr > 0)
            {
            	m_length = Double.valueOf(lenval);
            	lengthFound = true;
            }
		} 

		public void ProcessName()
		{
			int nchar = 0;
	        while ((m_newick[m_curr_char] != ':') &&
	        	   (m_newick[m_curr_char] != ',') &&
             	   (m_newick[m_curr_char] != ')'))
				{
	        		m_name += m_newick[m_curr_char];
					++m_curr_char;
					nchar++;
				}
	        if (nchar > maxNameLen)
	        {
	        	maxNameLen = nchar;
	        }
		} 	
	
		public void ProcessComment()
		// NB: We don't do anything with the contents of the comment yet
		{
	        while (m_newick[m_curr_char] != ']') 
			{
				++m_curr_char;
			}
		}
		
		public void DumpTree(NewickNode node)
		{
			// Debug dump of tree outward from a node to System.out
			System.out.println("Node number: " + node.m_nodenum);
			if (node.m_type == NodeType.ROOT)
			{
				System.out.println("Root node");
			}
			else if (node.m_type == NodeType.NODE)
			{
				System.out.println("Internal node");
			}
			else if (node.m_type == NodeType.TIP)
			{
				System.out.println("Tip node Name: " + node.m_name);
			}
			else
			{
				System.out.println("ERROR: Unknown node type");
			}
			System.out.println("Length: " + node.m_length);
			System.out.println("Xcor: " + node.m_corPt.x);
			System.out.println("Ycor: " + node.m_corPt.y);
			int n = 1;
			Iterator<NewickNode> itr = node.m_children.iterator();
			while(itr.hasNext())
			{
				System.out.println("     Child: "+ n );
				n++;
				DumpTree(itr.next());
			}
		}
		
		public void PlotTree(NewickNode node, Graphics g2d)
		{
			if(node.m_type == NodeType.TIP)
			{
				// print name
				node.m_dispPt.x = xOffset + rootLen + (int)(node.m_corPt.x * xScaleFactor);
				node.m_dispPt.y = (int)(node.m_corPt.y * yScaleFactor);
				int labelStartX = node.m_dispPt.x + nameOffsetX;
				int labelStartY = node.m_dispPt.y + nameTweakY;
				g2d.setColor(node.m_color);
				g2d.drawString(node.m_name,  labelStartX,  labelStartY);
				//System.out.println("Node Num: " + node.m_nodenum + " Type: " + node.m_type + " Xdisp: " + node.m_dispPt.x + " Ydisp: " + node.m_dispPt.y);
			
				return;
			}
			else if (node.m_type == NodeType.ROOT)
			{
				// draw root
				node.m_dispPt.x = xOffset + rootLen;
				node.m_dispPt.y = (int)(node.m_corPt.y * yScaleFactor);
				g2d.setColor(node.m_color);
				g2d.drawLine(xOffset, node.m_dispPt.y, node.m_dispPt.x,  node.m_dispPt.y); // root stub - arbitrary length
				
				// this is for the branch selector code
				newseg = new BranchSeg();
				newseg.m_startnode = node.m_parent;
				newseg.m_endnode = node;
				newseg.m_hbranch = new Line2D.Float(xOffset, node.m_dispPt.y, node.m_dispPt.x,  node.m_dispPt.y);
				branchSegs.add(newseg);				
				
			}
			
			Iterator<NewickNode> itr = node.m_children.iterator();
			while(itr.hasNext())
			{
				// draw branch
				NewickNode child = itr.next();
				child.m_dispPt.x = xOffset + rootLen + (int)(child.m_corPt.x * xScaleFactor);
				child.m_dispPt.y = (int)(child.m_corPt.y * yScaleFactor);
				g2d.setColor(child.m_color);
				g2d.drawLine(child.GetParent().m_dispPt.x, child.GetParent().m_dispPt.y, child.GetParent().m_dispPt.x,  child.m_dispPt.y); // y offset line
				g2d.drawLine(child.GetParent().m_dispPt.x, child.m_dispPt.y, child.m_dispPt.x,  child.m_dispPt.y); // x branch
				
				if (showLens)
				{
					Integer lenX = child.GetParent().m_dispPt.x + (child.m_dispPt.x - child.GetParent().m_dispPt.x)/2 - 10;
					g2d.drawString(Double.toString(child.m_length), lenX,  child.m_dispPt.y - 5);
				}

				g2d.setColor(Color.BLACK);
				//System.out.println("Node Num: " + node.m_nodenum + " Type: " + node.m_type + " Xdisp: " + node.m_dispPt.x + " Ydisp: " + node.m_dispPt.y);

				// this is for the branch selector code
				newseg = new BranchSeg();
				newseg.m_startnode = child.m_parent;
				newseg.m_endnode = child;
				newseg.m_vbranch = new Line2D.Float(child.GetParent().m_dispPt.x, child.GetParent().m_dispPt.y, child.GetParent().m_dispPt.x,  child.m_dispPt.y);
				newseg.m_hbranch = new Line2D.Float(child.GetParent().m_dispPt.x, child.m_dispPt.y, child.m_dispPt.x,  child.m_dispPt.y);
				branchSegs.add(newseg);	
				
				PlotTree(child, g2d);
			}
		}
		
		public void ClearTree(NewickNode node)
		{
			// reset color to black after edit
			if(node.m_type == NodeType.TIP)
			{
				node.m_color = Color.BLACK;
				return;
			}
			else if (node.m_type == NodeType.ROOT)
			{
				node.m_color = Color.BLACK;			
			}
			
			Iterator<NewickNode> itr = node.m_children.iterator();
			while(itr.hasNext())
			{
				NewickNode child = itr.next();
				node.m_color = Color.BLACK;
				ClearTree(child);
			}
		}
	}
	
	public class BranchSeg{
		Line2D.Float m_hbranch;
		Line2D.Float m_vbranch;
		NewickNode m_startnode;
		NewickNode m_endnode;
		public BranchSeg()
		{
			m_hbranch = null;
			m_vbranch = null;
			m_startnode = null;
			m_endnode = null;
		}
		
		public double distanceToBranch(Point2D.Float pick)
		{
			double a = pick.x - m_hbranch.x1;
			double b = pick.y - m_hbranch.y1;
			double c = m_hbranch.x2 - m_hbranch.x1;
			double d = m_hbranch.y2 - m_hbranch.y1;
			
			double dot_prod = a*c + b*d;
			double len_sq = c*c + d*d;
			double param = dot_prod / len_sq;
			
			double xx;
			double yy;
			
			if (param < 0 || (m_hbranch.x1 == m_hbranch.x2 && m_hbranch.y1 == m_hbranch.y2))
			{
				xx = m_hbranch.x1;
				yy = m_hbranch.y1;
			}
			else if (param > 1)
			{
				xx = m_hbranch.x2;
				yy = m_hbranch.y2;
			}
			else
			{
				xx = m_hbranch.x1 + param*c;
				yy = m_hbranch.y1 + param*d;
			}
			double dx = pick.x - xx;
			double dy = pick.y - yy;
			double dist = Math.sqrt(dx*dx + dy*dy);
			
			if (m_vbranch != null)
			{
				a = pick.x - m_vbranch.x1;
				b = pick.y - m_vbranch.y1;
				c = m_vbranch.x2 - m_vbranch.x1;
				d = m_vbranch.y2 - m_vbranch.y1;
				
				dot_prod = a*c + b*d;
				len_sq = c*c + d*d;
				param = dot_prod / len_sq;
				
				if (param < 0 || (m_vbranch.x1 == m_vbranch.x2 && m_vbranch.y1 == m_vbranch.y2))
				{
					xx = m_vbranch.x1;
					yy = m_vbranch.y1;
				}
				else if (param > 1)
				{
					xx = m_vbranch.x2;
					yy = m_vbranch.y2;
				}
				else
				{
					xx = m_vbranch.x1 + param*c;
					yy = m_vbranch.y1 + param*d;
				}
				dx = pick.x - xx;
				dy = pick.y - yy;
				if (Math.sqrt(dx*dx + dy*dy) < dist)
				{
					dist = Math.sqrt(dx*dx + dy*dy);
				}		
			}
			
			return dist;
		}
	}
}
