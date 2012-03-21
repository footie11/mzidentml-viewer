/*
 * MzIdentMLViewer.java
 *
 * Created on 15-Dec-2010, 10:34:37
 * 
 * Fawaz Ghali
 * 
 */
package mzidentmlviewer;

import com.compomics.util.gui.spectrum.DefaultSpectrumAnnotation;
import java.awt.GridLayout;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.awt.event.MouseEvent;
import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import java.util.regex.Matcher;
import java.util.regex.Pattern;
import javax.swing.JCheckBox;
import javax.swing.JFileChooser;
import javax.swing.JOptionPane;
import javax.swing.JScrollPane;
import javax.swing.ListSelectionModel;
import javax.swing.UIManager;
import javax.swing.UnsupportedLookAndFeelException;
import javax.swing.table.DefaultTableModel;
import org.jdesktop.swingx.JXTable;


import uk.ac.ebi.jmzidml.model.mzidml.CvParam;
import uk.ac.ebi.jmzidml.model.mzidml.DBSequence;
import uk.ac.ebi.jmzidml.model.mzidml.Enzyme;
import uk.ac.ebi.jmzidml.model.mzidml.Fragmentation;
import uk.ac.ebi.jmzidml.model.mzidml.IonType;
import uk.ac.ebi.jmzidml.model.mzidml.Modification;
import uk.ac.ebi.jmzidml.model.mzidml.Peptide;
import uk.ac.ebi.jmzidml.model.mzidml.PeptideEvidence;
import uk.ac.ebi.jmzidml.model.mzidml.ProteinAmbiguityGroup;
import uk.ac.ebi.jmzidml.model.mzidml.ProteinDetectionHypothesis;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentificationItem;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentificationList;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentificationResult;
import uk.ac.ebi.jmzidml.xml.io.MzIdentMLUnmarshaller;
import util.MzIdentMLFilter;
import util.ProgressBarDialog;
import com.compomics.util.gui.spectrum.SpectrumPanel;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Toolkit;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.zip.GZIPInputStream;
import javax.swing.text.SimpleAttributeSet;
import javax.swing.text.StyleConstants;
import javax.xml.bind.JAXBException;
import mzidentml.FalseDiscoveryRate;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import uk.ac.ebi.jmzidml.MzIdentMLElement;
import uk.ac.ebi.jmzidml.model.mzidml.AnalysisProtocolCollection;
import uk.ac.ebi.jmzidml.model.mzidml.AnalysisSoftware;
import uk.ac.ebi.jmzidml.model.mzidml.Enzymes;
import uk.ac.ebi.jmzidml.model.mzidml.ModificationParams;
import uk.ac.ebi.jmzidml.model.mzidml.Param;
import uk.ac.ebi.jmzidml.model.mzidml.PeptideEvidenceRef;
import uk.ac.ebi.jmzidml.model.mzidml.PeptideHypothesis;
import uk.ac.ebi.jmzidml.model.mzidml.ProteinDetectionProtocol;
import uk.ac.ebi.jmzidml.model.mzidml.SearchModification;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentificationProtocol;
import uk.ac.ebi.jmzidml.model.mzidml.Tolerance;
import util.CsvFileFilter;

import util.TheoreticalFragmentation;

/**
 *
 * @author Fawaz Ghali
 */
public class MzIdentMLViewer extends javax.swing.JFrame {
    // mzIdentML file

    private MzIdentMLFilter mzIdentMLFilter = new MzIdentMLFilter();
    private MzIdentMLUnmarshaller mzIdentMLUnmarshaller = null;
    // GUI tables
    private JXTable proteinAmbiguityGroupTable;
    private JXTable proteinDetectionHypothesisTable;
    private JXTable spectrumIdentificationItemProteinViewTable;
    private JXTable spectrumIdentificationResultTable;
    private JXTable spectrumIdentificationItemTable;
    private JXTable peptideEvidenceTable;
    private JXTable fragmentationTable;
    private JXTable dBSequenceTable;
    private JXTable protocolTable;
    //Hashmaps to store mzIdentML data
    private Fragmentation fragmentation;
    //GUI
    private ProgressBarDialog progressBarDialog;
    private JCheckBox[] filterCheckBoxIon = new JCheckBox[12];
    private JCheckBox[] filterCheckBoxCharge = new JCheckBox[12];
    private List<String> filterListIon = new ArrayList();
    private List<String> filterListCharge = new ArrayList();
    private List<IonType> ionTypeList = null;
    boolean newHeadersIIProteinTable = false;
    boolean newHeadersIRTable = false;
    private SpectrumPanel spectrumPanel;
    private Vector<DefaultSpectrumAnnotation> peakAnnotation = new Vector();
    //stats
    private FalseDiscoveryRate falseDiscoveryRate = null;
    private double falsePositiveSii;
    private double truePositiveSii;
    private double fdrSii;
    private double falsePositiveProteins;
    private double truePositiveProteins;
    private double fdrProteins;
    //Lists to store mzIdentML data
    private List<ProteinDetectionHypothesis> pDHListPassThreshold;
    private List<SpectrumIdentificationItem> spectrumIdentificationItemListForSpecificResult;
    // FDR and Stats
    private List<SpectrumIdentificationItem> sIIListPassThreshold = new ArrayList();
    private List<SpectrumIdentificationItem> sIIListBelowThreshold = new ArrayList();
    private List<SpectrumIdentificationItem> sIIListPassThresholdRankOne = new ArrayList();
    private List<Peptide> peptideListNonReduntant = new ArrayList();
    private List<SpectrumIdentificationItem> sIIListIsDecoyFalse = new ArrayList();
    private List<SpectrumIdentificationItem> sIIListIsDecoyTrue = new ArrayList();
    private boolean firstTab, secondTab, thirdTab, fourthTab, fifthTab;

    /**
     * Creates new form MzIdentMLViewer
     */
    public MzIdentMLViewer() {
        // Look and Feel
        try {
            UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
        } catch (ClassNotFoundException | InstantiationException | IllegalAccessException | UnsupportedLookAndFeelException ex) {
        }

        // Swing init components
        initComponents();

        // protein tab
        // protein Ambiguity Group Table
        proteinAmbiguityGroupTable = new JXTable() {
        };
        proteinAmbiguityGroupTable.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
        proteinAmbiguityGroupTable.addMouseListener(new java.awt.event.MouseAdapter() {

            @Override
            public void mouseClicked(java.awt.event.MouseEvent evt) {
                proteinAmbiguityGroupTableMouseClicked(evt);
            }
        });
        JScrollPane jProteinAmbiguityGroupScrollPane = new JScrollPane(proteinAmbiguityGroupTable);
        jProteinAmbiguityGroupPanel.setLayout(new java.awt.BorderLayout());
        jProteinAmbiguityGroupPanel.add(jProteinAmbiguityGroupScrollPane);
        proteinAmbiguityGroupTable.getTableHeader().setReorderingAllowed(false);
        // protein Detection Hypothesis Table
        proteinDetectionHypothesisTable = new JXTable() {
        };
        proteinDetectionHypothesisTable.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
        proteinDetectionHypothesisTable.addMouseListener(new java.awt.event.MouseAdapter() {

            @Override
            public void mouseClicked(java.awt.event.MouseEvent evt) {
                proteinDetectionHypothesisTableMouseClicked(evt);
            }
        });
        JScrollPane jProteinDetectionHypothesisPane = new JScrollPane(proteinDetectionHypothesisTable);
        jProteinDetectionHypothesisPanel.setLayout(new java.awt.BorderLayout());
        jProteinDetectionHypothesisPanel.add(jProteinDetectionHypothesisPane);
        proteinDetectionHypothesisTable.getTableHeader().setReorderingAllowed(false);
        // spectrum Identification Item Protein View Table
        spectrumIdentificationItemProteinViewTable = new JXTable() {
        };
        spectrumIdentificationItemProteinViewTable.setSortable(true);
        spectrumIdentificationItemProteinViewTable.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
        spectrumIdentificationItemProteinViewTable.addMouseListener(new java.awt.event.MouseAdapter() {

            @Override
            public void mouseClicked(java.awt.event.MouseEvent evt) {
                spectrumIdentificationItemProteinTableeMouseClicked(evt);
            }
        });
        JScrollPane jspectrumIdentificationItemProteinTableScrollPane = new JScrollPane(spectrumIdentificationItemProteinViewTable);
        jSpectrumIdentificationItemProteinPanel.setLayout(new java.awt.BorderLayout());
        jSpectrumIdentificationItemProteinPanel.add(jspectrumIdentificationItemProteinTableScrollPane);
        spectrumIdentificationItemProteinViewTable.getTableHeader().setReorderingAllowed(false);
        // spectrum tab
        // spectrum Identification Result Table
        spectrumIdentificationResultTable = new JXTable() {
        };
        spectrumIdentificationResultTable.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
        spectrumIdentificationResultTable.addMouseListener(new java.awt.event.MouseAdapter() {

            @Override
            public void mouseClicked(java.awt.event.MouseEvent evt) {
                spectrumIdentificationResultTableMouseClicked(evt);
            }
        });
        JScrollPane jSpectrumIdentificationResultTableScrollPane = new JScrollPane(spectrumIdentificationResultTable);
        jSpectrumIdentificationResultPanel.setLayout(new java.awt.BorderLayout());
        jSpectrumIdentificationResultPanel.add(jSpectrumIdentificationResultTableScrollPane);
        spectrumIdentificationResultTable.getTableHeader().setReorderingAllowed(false);
        //spectrum Identification Item Table
        spectrumIdentificationItemTable = new JXTable() {
        };

        JScrollPane jSpectrumIdentificationItemTableScrollPane = new JScrollPane(spectrumIdentificationItemTable);
        jSpectrumIdentificationItemPanel.setLayout(new java.awt.BorderLayout());
        jSpectrumIdentificationItemPanel.add(jSpectrumIdentificationItemTableScrollPane);
        spectrumIdentificationItemTable.getTableHeader().setReorderingAllowed(false);
        spectrumIdentificationItemTable.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
        spectrumIdentificationItemTable.addMouseListener(new java.awt.event.MouseAdapter() {

            @Override
            public void mouseClicked(java.awt.event.MouseEvent evt) {
                spectrumIdentificationItemTableMouseClicked(evt);
            }
        });

        //peptide Evidence Table 
        peptideEvidenceTable = new JXTable() {
        };

        JScrollPane jPeptideEvidenceTableScrollPane = new JScrollPane(peptideEvidenceTable);
        jPeptideEvidencePanel.setLayout(new java.awt.BorderLayout());
        jPeptideEvidencePanel.add(jPeptideEvidenceTableScrollPane);
        peptideEvidenceTable.getTableHeader().setReorderingAllowed(false);
        peptideEvidenceTable.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
        peptideEvidenceTable.addMouseListener(new java.awt.event.MouseAdapter() {

            @Override
            public void mouseClicked(java.awt.event.MouseEvent evt) {
                peptideEvidenceTableMouseClicked(evt);
            }
        });


        //fragmentation Table
        fragmentationTable = new JXTable() {
        };
        fragmentationTable.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
        fragmentationTable.addMouseListener(new java.awt.event.MouseAdapter() {

            @Override
            public void mouseClicked(java.awt.event.MouseEvent evt) {
                fragmentationTableMouseClicked(evt);
            }
        });
        JScrollPane jFragmentationTableScrollPane = new JScrollPane(fragmentationTable);
        jFragmentationPanel.setLayout(new java.awt.BorderLayout());
        jFragmentationPanel.add(jFragmentationTableScrollPane);
        fragmentationTable.getTableHeader().setReorderingAllowed(false);


        //dBSequenceTable Table
        dBSequenceTable = new JXTable() {
        };
        dBSequenceTable.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
        dBSequenceTable.addMouseListener(new java.awt.event.MouseAdapter() {

            @Override
            public void mouseClicked(java.awt.event.MouseEvent evt) {
                dBSequenceTableMouseClicked(evt);
            }
        });
        JScrollPane jdBSequenceTableScrollPane = new JScrollPane(dBSequenceTable);
        dBSequencePanel.setLayout(new java.awt.BorderLayout());
        dBSequencePanel.add(jdBSequenceTableScrollPane);
        dBSequenceTable.getTableHeader().setReorderingAllowed(false);


        //protocolTable Table
        protocolTable = new JXTable() {
        };
        protocolTable.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
        protocolTable.addMouseListener(new java.awt.event.MouseAdapter() {

            @Override
            public void mouseClicked(java.awt.event.MouseEvent evt) {
                protocolTableMouseClicked(evt);
            }
        });
        JScrollPane jprotocolTableScrollPane = new JScrollPane(protocolTable);
        protocolSummaryPanel.setLayout(new java.awt.BorderLayout());
        protocolSummaryPanel.add(jprotocolTableScrollPane);
        protocolTable.getTableHeader().setReorderingAllowed(false);
        // tables cannot be edited
        proteinAmbiguityGroupTable.setEditable(false);
        proteinDetectionHypothesisTable.setEditable(false);
        spectrumIdentificationItemProteinViewTable.setEditable(false);
        spectrumIdentificationResultTable.setEditable(false);
        spectrumIdentificationItemTable.setEditable(false);
        fragmentationTable.setEditable(false);
        peptideEvidenceTable.setEditable(false);
        dBSequenceTable.setEditable(false);
        spectrumIdentificationResultTable.setToolTipText("this corresponds to Spectrum Identification Result in mzIdentML");
        Dimension dim = Toolkit.getDefaultToolkit().getScreenSize();
        setSize(dim.width, dim.height - 40);
        setLocationRelativeTo(getRootPane());
        repaint();
    }

    private void protocolTableMouseClicked(MouseEvent evt) {
    }

    private void peptideEvidenceTableMouseClicked(MouseEvent evt) {
        int row = peptideEvidenceTable.getSelectedRow();
        if (row != -1) {
            String db_ref = (String) peptideEvidenceTable.getValueAt(row, 6);

            int rowCount = dBSequenceTable.getModel().getRowCount();
            for (int i = 0; i < rowCount; i++) {
                if (db_ref.equals((String) dBSequenceTable.getValueAt(i, 0))) {

                    dBSequenceTable.setRowSelectionInterval(i, i);
                }

            }
        }
        mainTabbedPane.setSelectedIndex(2);
    }

    private void dBSequenceTableMouseClicked(MouseEvent evt) {
    }

    public void createTables() {
        // dBSequence view dBSequenceTable
        String[] dBSequenceTableHeaders = new String[]{"ID", "Accession", "DB Ref", "Seq", "Protein Description"};
        dBSequenceTable.setModel(new javax.swing.table.DefaultTableModel(new Object[][]{}, dBSequenceTableHeaders) {
        });
        while (((DefaultTableModel) dBSequenceTable.getModel()).getRowCount() > 0) {
            ((DefaultTableModel) dBSequenceTable.getModel()).removeRow(0);
        }
        // protocolTable
        String[] protocolTableHeaders = new String[]{"", ""};
        protocolTable.setModel(new javax.swing.table.DefaultTableModel(new Object[][]{}, protocolTableHeaders) {
        });
        while (((DefaultTableModel) protocolTable.getModel()).getRowCount() > 0) {
            ((DefaultTableModel) protocolTable.getModel()).removeRow(0);
        }
        // protein view
        String[] proteinAmbiguityGroupTableHeaders = new String[]{"ID", "Name", "Protein Accessions"};
        String[] proteinDetectionHypothesisTableHeaders = new String[]{"ID", "Accession", "Scores", "P-values", "Number of peptides", "Is Decoy"};
        String[] spectrumIdentificationItemProteinViewTableHeaders = new String[]{"Peptide Sequence", "SII", "Name", "Score", "Expectation value"};
        proteinAmbiguityGroupTable.setModel(new javax.swing.table.DefaultTableModel(new Object[][]{}, proteinAmbiguityGroupTableHeaders) {
        });
        proteinDetectionHypothesisTable.setModel(new javax.swing.table.DefaultTableModel(new Object[][]{}, proteinDetectionHypothesisTableHeaders) {
        });
        spectrumIdentificationItemProteinViewTable.setModel(new javax.swing.table.DefaultTableModel(new Object[][]{}, spectrumIdentificationItemProteinViewTableHeaders) {
        });
        jProteinSequenceTextPane.setText("");
        jProteinSequenceTextPane.setText("");

        while (((DefaultTableModel) proteinAmbiguityGroupTable.getModel()).getRowCount() > 0) {
            ((DefaultTableModel) proteinAmbiguityGroupTable.getModel()).removeRow(0);
        }
        while (((DefaultTableModel) proteinDetectionHypothesisTable.getModel()).getRowCount() > 0) {
            ((DefaultTableModel) proteinDetectionHypothesisTable.getModel()).removeRow(0);
        }
        while (((DefaultTableModel) spectrumIdentificationItemProteinViewTable.getModel()).getRowCount() > 0) {
            ((DefaultTableModel) spectrumIdentificationItemProteinViewTable.getModel()).removeRow(0);
        }

        //spectrum view
        int spectrumIdentificationResultCvParamLengs = 0;
        int spectrumIdentificationItemCvParamLengs = 0;
        String[] sir = null;
        String[] sii = null;
        List<SpectrumIdentificationResult> spectrumIdentificationResultList;

        spectrumIdentificationResultList = new ArrayList();


        Iterator<SpectrumIdentificationList> iterspectrumIdentificationList = mzIdentMLUnmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.SpectrumIdentificationList);
        while (iterspectrumIdentificationList.hasNext()) {
            SpectrumIdentificationList spectrumIdentificationList = iterspectrumIdentificationList.next();


            for (int j = 0; j < spectrumIdentificationList.getSpectrumIdentificationResult().size(); j++) {
                spectrumIdentificationResultList.add(spectrumIdentificationList.getSpectrumIdentificationResult().get(j));
            }
        }

        if (spectrumIdentificationResultList.size() > 0) {

            for (int i = 0; i < spectrumIdentificationResultList.size(); i++) {
                SpectrumIdentificationResult spectrumIdentificationResult = spectrumIdentificationResultList.get(i);
                if (spectrumIdentificationResultCvParamLengs < spectrumIdentificationResult.getCvParam().size()) {
                    spectrumIdentificationResultCvParamLengs = spectrumIdentificationResult.getCvParam().size();
                    sir = new String[spectrumIdentificationResult.getCvParam().size()];
                    for (int j = 0; j < spectrumIdentificationResult.getCvParam().size(); j++) {
                        sir[j] = spectrumIdentificationResult.getCvParam().get(j).getName();
                    }
                }
                List<SpectrumIdentificationItem> spectrumIdentificationItemList = spectrumIdentificationResult.getSpectrumIdentificationItem();
                for (int j = 0; j < spectrumIdentificationItemList.size(); j++) {
                    SpectrumIdentificationItem spectrumIdentificationItem = spectrumIdentificationItemList.get(j);
                    if (spectrumIdentificationItemCvParamLengs < spectrumIdentificationItem.getCvParam().size()) {
                        spectrumIdentificationItemCvParamLengs = spectrumIdentificationItem.getCvParam().size();
                    }
                    sii = new String[spectrumIdentificationItem.getCvParam().size()];
                    for (int k = 0; k < spectrumIdentificationItem.getCvParam().size(); k++) {
                        sii[k] = spectrumIdentificationItem.getCvParam().get(k).getName();

                    }

                }
            }
        }

        String[] spectrumIdentificationResultTableHeaders = new String[spectrumIdentificationResultCvParamLengs + 3];
        spectrumIdentificationResultTableHeaders[0] = "ID";
        spectrumIdentificationResultTableHeaders[1] = "Spectrum ID";
        spectrumIdentificationResultTableHeaders[2] = "Spectra Data REF";
        if (sir != null) {
            for (int i = 0; i < sir.length; i++) {
                String string = sir[i];
                spectrumIdentificationResultTableHeaders[3 + i] = string;

            }
        }

        String[] spectrumIdentificationItemTableHeaders = new String[spectrumIdentificationItemCvParamLengs + 7];
        spectrumIdentificationItemTableHeaders[0] = "ID";
        spectrumIdentificationItemTableHeaders[1] = "Peptide Sequence";
        spectrumIdentificationItemTableHeaders[2] = "Modification";
        spectrumIdentificationItemTableHeaders[3] = "Calculated MassToCharge";
        spectrumIdentificationItemTableHeaders[4] = "Experimental MassToCharge";
        spectrumIdentificationItemTableHeaders[5] = "Rank";
        spectrumIdentificationItemTableHeaders[6] = "Is Decoy";
        if (sii != null) {
            for (int i = 0; i < sii.length; i++) {
                String string = sii[i];
                spectrumIdentificationItemTableHeaders[7 + i] = string;
            }
        }

//        String[] peptideEvidenceTableHeaders = new String[peptideEvidenceCvParamLengs + 7];
        String[] peptideEvidenceTableHeaders = new String[7];
        peptideEvidenceTableHeaders[0] = "Start";
        peptideEvidenceTableHeaders[1] = "End";
        peptideEvidenceTableHeaders[2] = "Pre";
        peptideEvidenceTableHeaders[3] = "Post";
        peptideEvidenceTableHeaders[4] = "IsDecoy";
        peptideEvidenceTableHeaders[5] = "Peptide Sequence";
        peptideEvidenceTableHeaders[6] = "dBSequence_ref";
//        if (pe != null) {
//            for (int i = 0; i < pe.length; i++) {
//                String string = pe[i];
//                peptideEvidenceTableHeaders[7 + i] = string;
//            }
//        }
//        

        String[] fragmentationTableHeaders = new String[]{"M/Z", "Intensity", "M Error", "Ion Type", "Charge"};

        spectrumIdentificationResultTable.setModel(new javax.swing.table.DefaultTableModel(new Object[][]{}, spectrumIdentificationResultTableHeaders) {
        });
        spectrumIdentificationItemTable.setModel(new javax.swing.table.DefaultTableModel(new Object[][]{}, spectrumIdentificationItemTableHeaders) {
        });
        peptideEvidenceTable.setModel(new javax.swing.table.DefaultTableModel(new Object[][]{}, peptideEvidenceTableHeaders) {
        });
        fragmentationTable.setModel(new javax.swing.table.DefaultTableModel(new Object[][]{}, fragmentationTableHeaders) {
        });
        while (((DefaultTableModel) spectrumIdentificationItemTable.getModel()).getRowCount() > 0) {
            ((DefaultTableModel) spectrumIdentificationItemTable.getModel()).removeRow(0);
        }
        while (((DefaultTableModel) spectrumIdentificationResultTable.getModel()).getRowCount() > 0) {
            ((DefaultTableModel) spectrumIdentificationResultTable.getModel()).removeRow(0);
        }
        while (((DefaultTableModel) peptideEvidenceTable.getModel()).getRowCount() > 0) {
            ((DefaultTableModel) peptideEvidenceTable.getModel()).removeRow(0);
        }
        while (((DefaultTableModel) fragmentationTable.getModel()).getRowCount() > 0) {
            ((DefaultTableModel) fragmentationTable.getModel()).removeRow(0);
        }



        //graph view
        while (jGraph.getComponents().length > 0) {
            jGraph.remove(0);
        }
        while (jExperimentalFilterPanel.getComponents().length > 0) {
            jExperimentalFilterPanel.remove(0);
        }

        jProteinDescriptionEditorPane.setText("");
        jGraph.validate();
        jGraph.repaint();

    }

    private void loadProteinAmbiguityGroupTable() {
        this.setCursor(new java.awt.Cursor(java.awt.Cursor.WAIT_CURSOR));
        String protein_accessions = "";
        pDHListPassThreshold = new ArrayList();

        Iterator<ProteinAmbiguityGroup> iterProteinAmbiguityGroup = mzIdentMLUnmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.ProteinAmbiguityGroup);
        List<ProteinDetectionHypothesis> proteinDetectionHypothesisList;
        while (iterProteinAmbiguityGroup.hasNext()) {
            ProteinAmbiguityGroup proteinAmbiguityGroup = iterProteinAmbiguityGroup.next();


            proteinDetectionHypothesisList = proteinAmbiguityGroup.getProteinDetectionHypothesis();

            if (proteinDetectionHypothesisList.size() > 0) {
                for (int j = 0; j < proteinDetectionHypothesisList.size(); j++) {
                    ProteinDetectionHypothesis proteinDetectionHypothesis = proteinDetectionHypothesisList.get(j);
                    if (proteinDetectionHypothesis.getId().startsWith("PDH_")) {
                        protein_accessions = protein_accessions + proteinDetectionHypothesis.getId().substring(4) + ";";
                    }

                    if (proteinDetectionHypothesis.isPassThreshold()) {
                        pDHListPassThreshold.add(proteinDetectionHypothesis);
                    }

                }
            }
            ((DefaultTableModel) proteinAmbiguityGroupTable.getModel()).addRow(new String[]{
                        proteinAmbiguityGroup.getId(),
                        proteinAmbiguityGroup.getName(),
                        protein_accessions
                    });
        }

        this.setCursor(new java.awt.Cursor(java.awt.Cursor.DEFAULT_CURSOR));
        proteinAmbiguityGroupTable.setSortable(true);
    }

    private void loadSpectrumIdentificationResultTable() {
        this.setCursor(new java.awt.Cursor(java.awt.Cursor.WAIT_CURSOR));
        Iterator<SpectrumIdentificationResult> iterSpectrumIdentificationResult = mzIdentMLUnmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.SpectrumIdentificationResult);
        while (iterSpectrumIdentificationResult.hasNext()) {


            SpectrumIdentificationResult spectrumIdentificationResult = iterSpectrumIdentificationResult.next();
            ((DefaultTableModel) spectrumIdentificationResultTable.getModel()).addRow(new String[]{
                        spectrumIdentificationResult.getId(),
                        spectrumIdentificationResult.getSpectrumID(),
                        spectrumIdentificationResult.getSpectraDataRef()
                    });
            List<CvParam> cvParamListspectrumIdentificationResult = spectrumIdentificationResult.getCvParam();

            for (int s = 0; s < cvParamListspectrumIdentificationResult.size(); s++) {
                CvParam cvParam = cvParamListspectrumIdentificationResult.get(s);
                String name = cvParam.getName();
                for (int j = 0; j < spectrumIdentificationResultTable.getModel().getColumnCount(); j++) {
                    if (spectrumIdentificationResultTable.getModel().getColumnName(j).equals(name)) {
                        ((DefaultTableModel) spectrumIdentificationResultTable.getModel()).setValueAt(cvParam.getValue(), ((DefaultTableModel) spectrumIdentificationResultTable.getModel()).getRowCount() - 1, j);
                    }
                }

            }

        }

        spectrumIdentificationResultTable.setSortable(true);

        this.setCursor(new java.awt.Cursor(java.awt.Cursor.DEFAULT_CURSOR));

    }

    /**
     * Creates fragmentation Table Mouse Clicked
     */
    private void fragmentationTableMouseClicked(MouseEvent evt) {
    }

    /**
     * Creates spectrum Identification Result Table Mouse Clicked
     */
    private void spectrumIdentificationResultTableMouseClicked(MouseEvent evt) {

        spectrumIdentificationItemTable.setSortable(false);
        this.setCursor(new java.awt.Cursor(java.awt.Cursor.WAIT_CURSOR));
        int row = spectrumIdentificationResultTable.getSelectedRow();
        if (row != -1) {
            try {
                while (((DefaultTableModel) spectrumIdentificationItemTable.getModel()).getRowCount() > 0) {
                    ((DefaultTableModel) spectrumIdentificationItemTable.getModel()).removeRow(0);
                }
                while (((DefaultTableModel) peptideEvidenceTable.getModel()).getRowCount() > 0) {
                    ((DefaultTableModel) peptideEvidenceTable.getModel()).removeRow(0);
                }
                while (((DefaultTableModel) fragmentationTable.getModel()).getRowCount() > 0) {
                    ((DefaultTableModel) fragmentationTable.getModel()).removeRow(0);
                }
                while (jGraph.getComponents().length > 0) {
                    jGraph.remove(0);
                }
                jGraph.validate();
                jGraph.repaint();
                spectrumIdentificationItemTable.scrollRowToVisible(0);
                String sir_id = (String) spectrumIdentificationResultTable.getModel().getValueAt(row, 0);
                SpectrumIdentificationResult spectrumIdentificationResult = mzIdentMLUnmarshaller.unmarshal(SpectrumIdentificationResult.class, sir_id);
                spectrumIdentificationItemListForSpecificResult = spectrumIdentificationResult.getSpectrumIdentificationItem();
                if (spectrumIdentificationItemListForSpecificResult.size() > 0) {

                    for (int i = 0; i < spectrumIdentificationItemListForSpecificResult.size(); i++) {
                        try {
                            SpectrumIdentificationItem spectrumIdentificationItem = spectrumIdentificationItemListForSpecificResult.get(i);
                            boolean isDecoy = checkIfSpectrumIdentificationItemIsDecoy(spectrumIdentificationItem);


                            Peptide peptide = mzIdentMLUnmarshaller.unmarshal(Peptide.class, spectrumIdentificationItem.getPeptideRef());
                            if (peptide != null) {
                                List<Modification> modificationList = peptide.getModification();
                                Modification modification;
                                String residues = null;
                                Integer location = null;
                                String modificationName = null;
                                CvParam modificationCvParam;
                                String combine = null;
                                if (modificationList.size() > 0) {
                                    modification = modificationList.get(0);
                                    location = modification.getLocation();
                                    if (modification.getResidues().size() > 0) {
                                        residues = modification.getResidues().get(0);
                                    }
                                    List<CvParam> modificationCvParamList = modification.getCvParam();
                                    if (modificationCvParamList.size() > 0) {
                                        modificationCvParam = modificationCvParamList.get(0);
                                        modificationName = modificationCvParam.getName();
                                    }
                                }
                                if (modificationName != null) {
                                    combine = modificationName;
                                }
                                if (residues != null) {
                                    combine = combine + " on residues: " + residues;
                                }
                                if (location != null) {
                                    combine = combine + " at location: " + location;
                                }
                                double calculatedMassToCharge = 0;
                                if (spectrumIdentificationItem.getCalculatedMassToCharge() != null) {
                                    calculatedMassToCharge = spectrumIdentificationItem.getCalculatedMassToCharge().doubleValue();
                                }
                                ((DefaultTableModel) spectrumIdentificationItemTable.getModel()).addRow(new Object[]{
                                            spectrumIdentificationItem.getId(),
                                            peptide.getPeptideSequence(),
                                            combine,
                                            roundTwoDecimals(calculatedMassToCharge),
                                            roundTwoDecimals(spectrumIdentificationItem.getExperimentalMassToCharge()),
                                            Integer.valueOf(spectrumIdentificationItem.getRank()),
                                            isDecoy
                                        });



                                List<CvParam> cvParamListspectrumIdentificationItem = spectrumIdentificationItem.getCvParam();

                                for (int s = 0; s < cvParamListspectrumIdentificationItem.size(); s++) {
                                    CvParam cvParam = cvParamListspectrumIdentificationItem.get(s);
                                    String accession = cvParam.getAccession();

                                    if (cvParam.getName().equals("peptide unique to one protein")) {
                                        ((DefaultTableModel) spectrumIdentificationItemTable.getModel()).setValueAt(1, ((DefaultTableModel) spectrumIdentificationItemTable.getModel()).getRowCount() - 1, 7 + s);
                                    } else if (accession.equals("MS:1001330")
                                            || accession.equals("MS:1001172")
                                            || accession.equals("MS:1001159")
                                            || accession.equals("MS:1001328")) {
                                        ((DefaultTableModel) spectrumIdentificationItemTable.getModel()).setValueAt(roundScientificNumbers(Double.valueOf(cvParam.getValue()).doubleValue()), ((DefaultTableModel) spectrumIdentificationItemTable.getModel()).getRowCount() - 1, 7 + s);
                                    } else {
                                        ((DefaultTableModel) spectrumIdentificationItemTable.getModel()).setValueAt(cvParam.getValue(), ((DefaultTableModel) spectrumIdentificationItemTable.getModel()).getRowCount() - 1, 7 + s);
                                    }
                                }

                            }
                        } catch (JAXBException ex) {
                            Logger.getLogger(MzIdentMLViewer.class.getName()).log(Level.SEVERE, null, ex);
                        }
                    }
                }
            } catch (JAXBException ex) {
                Logger.getLogger(MzIdentMLViewer.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
        spectrumIdentificationItemTable.setSortable(true);

        this.setCursor(new java.awt.Cursor(java.awt.Cursor.DEFAULT_CURSOR));
    }

    private double roundTwoDecimals(double d) {
        DecimalFormat twoDForm = new DecimalFormat("#.##");
        return Double.valueOf(twoDForm.format(d));
    }

    private double roundThreeDecimals(double d) {
        DecimalFormat twoDForm = new DecimalFormat("#.###");
        return Double.valueOf(twoDForm.format(d));
    }

    private String roundScientificNumbers(double d) {
        DecimalFormat formatter = new DecimalFormat("0.#E00");
        String f = formatter.format(d);
        return f;
    }

    private boolean checkIfSpectrumIdentificationItemIsDecoy(SpectrumIdentificationItem spectrumIdentificationItem) {
        boolean result = false;

        List<PeptideEvidenceRef> peptideEvidenceRefList = spectrumIdentificationItem.getPeptideEvidenceRef();
        for (int i = 0; i < peptideEvidenceRefList.size(); i++) {
            try {
                PeptideEvidenceRef peptideEvidenceRef = peptideEvidenceRefList.get(i);
                PeptideEvidence peptiedEvidence = mzIdentMLUnmarshaller.unmarshal(PeptideEvidence.class, peptideEvidenceRef.getPeptideEvidenceRef());
                if (peptiedEvidence != null && peptiedEvidence.isIsDecoy()) {
                    result = true;
                    break;
                }
            } catch (JAXBException ex) {
                Logger.getLogger(MzIdentMLViewer.class.getName()).log(Level.SEVERE, null, ex);
            }


        }
        return result;
    }

    /**
     * Creates spectrum Identification Item Table Mouse Clicked
     */
    private void spectrumIdentificationItemTableMouseClicked(MouseEvent evt) {
        fragmentationTable.setSortable(false);

        this.setCursor(new java.awt.Cursor(java.awt.Cursor.WAIT_CURSOR));
        int row = spectrumIdentificationItemTable.getSelectedRow();
        if (row != -1) {
            while (((DefaultTableModel) fragmentationTable.getModel()).getRowCount() > 0) {
                ((DefaultTableModel) fragmentationTable.getModel()).removeRow(0);
            }
            while (((DefaultTableModel) peptideEvidenceTable.getModel()).getRowCount() > 0) {
                ((DefaultTableModel) peptideEvidenceTable.getModel()).removeRow(0);
            }

            fragmentationTable.scrollRowToVisible(0);
            SpectrumIdentificationItem spectrumIdentificationItem = spectrumIdentificationItemListForSpecificResult.get(row);

            if (spectrumIdentificationItem != null) {
                List<PeptideEvidenceRef> peptideEvidenceRefList = spectrumIdentificationItem.getPeptideEvidenceRef();
                if (peptideEvidenceRefList != null) {
                    for (int i = 0; i < peptideEvidenceRefList.size(); i++) {
                        try {
                            PeptideEvidenceRef peptideEvidenceRef = peptideEvidenceRefList.get(i);

                            PeptideEvidence peptideEvidence = mzIdentMLUnmarshaller.unmarshal(PeptideEvidence.class, peptideEvidenceRef.getPeptideEvidenceRef());

                            ((DefaultTableModel) peptideEvidenceTable.getModel()).addRow(new Object[]{
                                        peptideEvidence.getStart(),
                                        peptideEvidence.getEnd(),
                                        peptideEvidence.getPre(),
                                        peptideEvidence.getPost(),
                                        peptideEvidence.isIsDecoy(),
                                        peptideEvidence.getPeptideRef(),
                                        peptideEvidence.getDBSequenceRef()
                                    // "<html><a href=>" +peptideEvidence.getDBSequenceRef()+"</a>"
                                    });
                        } catch (JAXBException ex) {
                            Logger.getLogger(MzIdentMLViewer.class.getName()).log(Level.SEVERE, null, ex);
                        }
                    }
                }
            }

            fragmentation = spectrumIdentificationItem.getFragmentation();


            if (fragmentation != null) {
                ionTypeList = fragmentation.getIonType();
                if (ionTypeList != null) {
                    for (int i = 0; i < ionTypeList.size(); i++) {
                        IonType ionType = ionTypeList.get(i);
                        CvParam cvParam = ionType.getCvParam();
                        if (!filterListIon.contains(cvParam.getName())) {
                            filterListIon.add(cvParam.getName());
                        }
                        if (!filterListCharge.contains(String.valueOf(ionType.getCharge()))) {
                            filterListCharge.add(String.valueOf(ionType.getCharge()));
                        }
                        List m_mz, m_intensity, m_error;
                        m_mz = ionType.getFragmentArray().get(0).getValues();
                        m_intensity = ionType.getFragmentArray().get(1).getValues();
                        m_error = ionType.getFragmentArray().get(2).getValues();
                        String type = cvParam.getName();
                        if (m_mz != null && !m_mz.isEmpty()) {
                            for (int j = 0; j < m_mz.size(); j++) {
                                ((DefaultTableModel) fragmentationTable.getModel()).addRow(new Object[]{
                                            Double.valueOf(m_mz.get(j).toString()),
                                            Double.valueOf(m_intensity.get(j).toString()),
                                            Double.valueOf(m_error.get(j).toString()),
                                            type,
                                            Integer.valueOf(ionType.getCharge())
                                        });


                            }
                        }
                    }
                }
                while (jExperimentalFilterPanel.getComponents().length > 0) {
                    jExperimentalFilterPanel.remove(0);
                }
                jExperimentalFilterPanel.validate();
                jExperimentalFilterPanel.repaint();
                if (filterListIon.isEmpty() && filterListCharge.isEmpty()) {
                    jExperimentalFilterPanel.setLayout(new GridLayout(0, 1));
                } else {
                    jExperimentalFilterPanel.setLayout(new GridLayout(0, filterListIon.size() + filterListCharge.size()));
                }

                for (int k = 0; k < filterListIon.size(); k++) {
                    String name = filterListIon.get(k);
                    filterCheckBoxIon[k] = new JCheckBox(name);
                    filterCheckBoxIon[k].setSelected(true);
                    filterCheckBoxIon[k].addItemListener(new ItemListener() {

                        @Override
                        public void itemStateChanged(ItemEvent E) {
                            updateGraph();
                        }
                    });
                    String type = name;
                    type = type.replaceFirst("frag:", "");
                    type = type.replaceFirst("ion", "");
                    type = type.replaceFirst("internal", "");

                    filterCheckBoxIon[k].setText(type);
                    jExperimentalFilterPanel.add(filterCheckBoxIon[k]);
                }

                for (int k = 0; k < filterListCharge.size(); k++) {
                    String name = filterListCharge.get(k);
                    filterCheckBoxCharge[k] = new JCheckBox(name);
                    filterCheckBoxCharge[k].setSelected(true);
                    filterCheckBoxCharge[k].addItemListener(new ItemListener() {

                        @Override
                        public void itemStateChanged(ItemEvent E) {
                            updateGraph();
                        }
                    });
                    String type = name;


                    filterCheckBoxCharge[k].setText(type);
                    jExperimentalFilterPanel.add(filterCheckBoxCharge[k]);
                }

                jExperimentalFilterPanel.repaint();
                jExperimentalFilterPanel.revalidate();

                double[] mzValuesAsDouble = new double[fragmentationTable.getModel().getRowCount()];
                double[] intensityValuesAsDouble = new double[fragmentationTable.getModel().getRowCount()];
                double[] m_errorValuesAsDouble = new double[fragmentationTable.getModel().getRowCount()];
                peakAnnotation.clear();
                for (int k = 0; k < fragmentationTable.getModel().getRowCount(); k++) {
                    mzValuesAsDouble[k] = (Double) (fragmentationTable.getModel().getValueAt(k, 0));

                    intensityValuesAsDouble[k] = (Double) (fragmentationTable.getModel().getValueAt(k, 1));
                    m_errorValuesAsDouble[k] = (Double) (fragmentationTable.getModel().getValueAt(k, 2));

                    String type = (String) fragmentationTable.getModel().getValueAt(k, 3);
                    type = type.replaceFirst("frag:", "");
                    type = type.replaceFirst("ion", "");
                    type = type.replaceFirst("internal", "");

                    peakAnnotation.add(
                            new DefaultSpectrumAnnotation(
                            mzValuesAsDouble[k],
                            m_errorValuesAsDouble[k],
                            Color.blue,
                            type));
                }
                if (mzValuesAsDouble.length > 0) {
                    spectrumPanel = new SpectrumPanel(
                            mzValuesAsDouble,
                            intensityValuesAsDouble,
                            spectrumIdentificationItem.getExperimentalMassToCharge(),
                            String.valueOf(spectrumIdentificationItem.getChargeState()),
                            spectrumIdentificationItem.getName());

                    spectrumPanel.setAnnotations(peakAnnotation);


                    while (jGraph.getComponents().length > 0) {
                        jGraph.remove(0);
                    }
                    jGraph.setLayout(new java.awt.BorderLayout());
                    jGraph.setLayout(new javax.swing.BoxLayout(jGraph, javax.swing.BoxLayout.LINE_AXIS));
                    jGraph.add(spectrumPanel);
                    jGraph.validate();
                    jGraph.repaint();
                    this.repaint();
                }
            }
        }
        fragmentationTable.setSortable(true);
        this.setCursor(new java.awt.Cursor(java.awt.Cursor.DEFAULT_CURSOR));
    }

    /**
     * Update Graph
     */
    private void updateGraph() {
        setCursor(new java.awt.Cursor(java.awt.Cursor.WAIT_CURSOR));
        int row = spectrumIdentificationItemTable.getSelectedRow();
        if (row != -1) {
            try {
                SpectrumIdentificationItem spectrumIdentificationItem = spectrumIdentificationItemListForSpecificResult.get(row);
                while (((DefaultTableModel) fragmentationTable.getModel()).getRowCount() > 0) {
                    ((DefaultTableModel) fragmentationTable.getModel()).removeRow(0);
                }
                for (int k = 0; k < filterListIon.size(); k++) {
                    String nameIon = filterListIon.get(k);
                    boolean isSelectedIon = filterCheckBoxIon[k].isSelected();
                    for (int z = 0; z < filterListCharge.size(); z++) {
                        String nameCharge = filterListCharge.get(z);
                        boolean isSelectedCharge = filterCheckBoxCharge[z].isSelected();

                        for (int i = 0; i < ionTypeList.size(); i++) {
                            IonType ionType = ionTypeList.get(i);
                            CvParam cvParam = ionType.getCvParam();
                            if ((nameIon.equals(cvParam.getName()) && isSelectedIon) && (nameCharge.equals(String.valueOf(ionType.getCharge())) && isSelectedCharge)) {
                                List m_mz, m_intensity, m_error;
                                m_mz = ionType.getFragmentArray().get(0).getValues();
                                m_intensity = ionType.getFragmentArray().get(1).getValues();
                                m_error = ionType.getFragmentArray().get(2).getValues();
                                if (m_mz != null && !m_mz.isEmpty()) {
                                    for (int j = 0; j < m_mz.size(); j++) {
                                        String type = cvParam.getName();
                                        ((DefaultTableModel) fragmentationTable.getModel()).addRow(new Object[]{
                                                    Double.valueOf(m_mz.get(j).toString()),
                                                    Double.valueOf(m_intensity.get(j).toString()),
                                                    Double.valueOf(m_error.get(j).toString()),
                                                    type,
                                                    Integer.valueOf(ionType.getCharge())
                                                });
                                    }
                                }
                            }
                        }
                    }
                }
                // Start of theortical values
                Peptide peptide = mzIdentMLUnmarshaller.unmarshal(Peptide.class, spectrumIdentificationItem.getPeptideRef());
                String find = peptide.getPeptideSequence();

                TheoreticalFragmentation tf = new TheoreticalFragmentation(find);

                List<Double> tmp = new ArrayList();
                if (bCheckBox.isSelected() && oneCheckBox.isSelected()) {
                    tmp.clear();
                    tmp = tf.getBIons(find, "1");
                    for (int j = 0; j < tmp.size(); j++) {

                        ((DefaultTableModel) fragmentationTable.getModel()).addRow(new Object[]{
                                    Double.valueOf(tmp.get(j).toString()),
                                    Double.valueOf("100"),
                                    Double.valueOf("0.0"),
                                    "T b+",
                                    +1
                                });
                    }


                }
                if (bCheckBox.isSelected() && twoCheckBox.isSelected()) {

                    tmp.clear();
                    tmp = tf.getBIons(find, "2");
                    for (int j = 0; j < tmp.size(); j++) {

                        ((DefaultTableModel) fragmentationTable.getModel()).addRow(new Object[]{
                                    Double.valueOf(tmp.get(j).toString()),
                                    Double.valueOf("100"),
                                    Double.valueOf("0.0"),
                                    "T b++",
                                    +1
                                });
                    }
                }
                if (yCheckBox.isSelected() && oneCheckBox.isSelected()) {


                    tmp.clear();
                    tmp = tf.getYIons(find, "1");
                    for (int j = 0; j < tmp.size(); j++) {

                        ((DefaultTableModel) fragmentationTable.getModel()).addRow(new Object[]{
                                    Double.valueOf(tmp.get(j).toString()),
                                    Double.valueOf("100"),
                                    Double.valueOf("0.0"),
                                    "T y+",
                                    +1
                                });
                    }
                }
                if (yCheckBox.isSelected() && twoCheckBox.isSelected()) {


                    tmp.clear();
                    tmp = tf.getYIons(find, "2");
                    for (int j = 0; j < tmp.size(); j++) {

                        ((DefaultTableModel) fragmentationTable.getModel()).addRow(new Object[]{
                                    Double.valueOf(tmp.get(j).toString()),
                                    Double.valueOf("100"),
                                    Double.valueOf("0.0"),
                                    "T y++",
                                    +1
                                });
                    }
                }



                // End of theortical values


                peakAnnotation.clear();
                double[] mzValuesAsDouble = new double[fragmentationTable.getModel().getRowCount()];
                double[] intensityValuesAsDouble = new double[fragmentationTable.getModel().getRowCount()];
                double[] m_errorValuesAsDouble = new double[fragmentationTable.getModel().getRowCount()];
                peakAnnotation.clear();
                for (int k = 0; k < fragmentationTable.getModel().getRowCount(); k++) {
                    mzValuesAsDouble[k] = (Double) (fragmentationTable.getModel().getValueAt(k, 0));
                    intensityValuesAsDouble[k] = (Double) (fragmentationTable.getModel().getValueAt(k, 1));
                    m_errorValuesAsDouble[k] = (Double) (fragmentationTable.getModel().getValueAt(k, 2));
                    String type = (String) fragmentationTable.getModel().getValueAt(k, 3);
                    type = type.replaceFirst("frag:", "");
                    type = type.replaceFirst("ion", "");
                    type = type.replaceFirst("internal", "");
                    peakAnnotation.add(
                            
                            
                            new DefaultSpectrumAnnotation(
                            mzValuesAsDouble[k],
                            m_errorValuesAsDouble[k],
                            Color.blue,
                            type));                                    // the annotation label
                }
                if (fragmentationTable.getModel().getRowCount() > 0) {
                    spectrumPanel = new SpectrumPanel(
                            mzValuesAsDouble,
                            intensityValuesAsDouble,
                            spectrumIdentificationItem.getExperimentalMassToCharge(),
                            String.valueOf(spectrumIdentificationItem.getChargeState()),
                            spectrumIdentificationItem.getName());
                    spectrumPanel.setAnnotations(peakAnnotation);
                    while (jGraph.getComponents().length > 0) {
                        jGraph.remove(0);
                    }
                    jGraph.setLayout(new java.awt.BorderLayout());
                    jGraph.setLayout(new javax.swing.BoxLayout(jGraph, javax.swing.BoxLayout.LINE_AXIS));
                    jGraph.add(spectrumPanel);
                    jGraph.validate();
                    jGraph.repaint();
                    this.repaint();
                }
                setCursor(new java.awt.Cursor(java.awt.Cursor.DEFAULT_CURSOR));
            } catch (JAXBException ex) {
                Logger.getLogger(MzIdentMLViewer.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
    }

    /**
     * Creates spectrum Identification Item Protein Table Mouse Clicked
     */
    private void spectrumIdentificationItemProteinTableeMouseClicked(MouseEvent evt) {
        this.setCursor(new java.awt.Cursor(java.awt.Cursor.WAIT_CURSOR));

        String sii_ref = (String) spectrumIdentificationItemProteinViewTable.getValueAt(spectrumIdentificationItemProteinViewTable.getSelectedRow(), 1);

        for (int i = 0; i < spectrumIdentificationResultTable.getRowCount(); i++) {
            try {
                String sir_id = (String) spectrumIdentificationResultTable.getValueAt(i, 0);

                SpectrumIdentificationResult sir = mzIdentMLUnmarshaller.unmarshal(SpectrumIdentificationResult.class, sir_id);
                List<SpectrumIdentificationItem> siiList = sir.getSpectrumIdentificationItem();
                for (int j = 0; j < siiList.size(); j++) {
                    SpectrumIdentificationItem spectrumIdentificationItem = siiList.get(j);
                    if (sii_ref.equals(spectrumIdentificationItem.getId())) {

                        spectrumIdentificationResultTable.setRowSelectionInterval(i, i);
                        spectrumIdentificationResultTableMouseClicked(evt);

                        spectrumIdentificationItemTable.setRowSelectionInterval(j, j);
                        spectrumIdentificationItemTableMouseClicked(evt);
                        break;

                    }

                }
            } catch (JAXBException ex) {
                Logger.getLogger(MzIdentMLViewer.class.getName()).log(Level.SEVERE, null, ex);
            }


        }

        this.setCursor(new java.awt.Cursor(java.awt.Cursor.DEFAULT_CURSOR));
        mainTabbedPane.setSelectedIndex(1);
    }

    /**
     * This method is called from within the constructor to initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is always
     * regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        fileChooser = new javax.swing.JFileChooser();
        mainPanel = new javax.swing.JPanel();
        mainTabbedPane = new javax.swing.JTabbedPane();
        proteinViewPanel = new javax.swing.JPanel();
        jProteinAmbiguityGroupPanel = new javax.swing.JPanel();
        jProteinDetectionHypothesisPanel = new javax.swing.JPanel();
        jProteinInfoPanel = new javax.swing.JPanel();
        jProteinSequencePanel = new javax.swing.JPanel();
        jProteinSequenceScrollPane = new javax.swing.JScrollPane();
        jProteinSequenceTextPane = new javax.swing.JTextPane();
        jSpectrumIdentificationItemProteinPanel = new javax.swing.JPanel();
        jScientificNameLabel = new javax.swing.JLabel();
        jScientificNameValueLabel = new javax.swing.JLabel();
        jProteinDescriptionPanel = new javax.swing.JPanel();
        jProteinDescriptionScrollPane = new javax.swing.JScrollPane();
        jProteinDescriptionEditorPane = new javax.swing.JEditorPane();
        spectrumViewPanel = new javax.swing.JPanel();
        jSpectrumIdentificationResultPanel = new javax.swing.JPanel();
        jSpectrumIdentificationItemPanel = new javax.swing.JPanel();
        jPeptideEvidencePanel = new javax.swing.JPanel();
        jSpectrumPanel = new javax.swing.JPanel();
        jFragmentationPanel = new javax.swing.JPanel();
        jGraph = new javax.swing.JPanel();
        jExperimentalFilterPanel = new javax.swing.JPanel();
        jTheoreticalFilterPanel = new javax.swing.JPanel();
        bCheckBox = new javax.swing.JCheckBox();
        yCheckBox = new javax.swing.JCheckBox();
        oneCheckBox = new javax.swing.JCheckBox();
        twoCheckBox = new javax.swing.JCheckBox();
        proteinDBViewPanel = new javax.swing.JPanel();
        dBSequencePanel = new javax.swing.JPanel();
        globalStatisticsPanel = new javax.swing.JPanel();
        summaryPanel = new javax.swing.JPanel();
        totalSIRLabel = new javax.swing.JLabel();
        totalSIRLabelValue = new javax.swing.JLabel();
        totalSIILabel = new javax.swing.JLabel();
        totalSIIbelowThresholdLabel = new javax.swing.JLabel();
        totalSIIbelowThresholdLabelValue = new javax.swing.JLabel();
        totalSIIaboveThresholdLabel = new javax.swing.JLabel();
        totalSIIaboveThresholdLabelValue = new javax.swing.JLabel();
        totalSIIaboveThresholdRankOneLabel = new javax.swing.JLabel();
        totalSIIaboveThresholdRankOneLabelValue = new javax.swing.JLabel();
        percentIdentifiedSpectraLabel = new javax.swing.JLabel();
        percentIdentifiedSpectraLabelValue = new javax.swing.JLabel();
        totalPeptidesaboveThresholdLabel = new javax.swing.JLabel();
        totalPeptidesaboveThresholdLabelValue = new javax.swing.JLabel();
        totalPAGsLabel = new javax.swing.JLabel();
        totalPAGsLabelValue = new javax.swing.JLabel();
        totalPDHsLabel = new javax.swing.JLabel();
        totalPDHsaboveThresholdLabel = new javax.swing.JLabel();
        totalPDHsaboveThresholdLabelValue = new javax.swing.JLabel();
        jSeparator1 = new javax.swing.JSeparator();
        isDecoySii = new javax.swing.JLabel();
        isDecoySiiValue = new javax.swing.JLabel();
        isDecoySiiFalse = new javax.swing.JLabel();
        isDecoySiiFalseValue = new javax.swing.JLabel();
        manualDecoyLabel = new javax.swing.JLabel();
        manualDecoyPrefix = new javax.swing.JLabel();
        manualDecoyRatio = new javax.swing.JLabel();
        manualDecoyPrefixValue = new javax.swing.JTextField();
        manualDecoyRatioValue = new javax.swing.JTextField();
        manualDecoy = new javax.swing.JCheckBox();
        fpSiiLabel = new javax.swing.JLabel();
        fpSiiValue = new javax.swing.JLabel();
        tpSiiLabel = new javax.swing.JLabel();
        tpSiiValue = new javax.swing.JLabel();
        fdrSiiLabel = new javax.swing.JLabel();
        fdrSiiValue = new javax.swing.JLabel();
        fpProteinsLabel = new javax.swing.JLabel();
        tpProteinsLabel = new javax.swing.JLabel();
        fdrProteinsLabel = new javax.swing.JLabel();
        fpProteinsValue = new javax.swing.JLabel();
        tpProteinsValue = new javax.swing.JLabel();
        fdrProteinsValue = new javax.swing.JLabel();
        manualCalculate = new javax.swing.JButton();
        jSeparator2 = new javax.swing.JSeparator();
        fdrPanel = new javax.swing.JPanel();
        thresholdLabel = new javax.swing.JLabel();
        thresholdValue = new javax.swing.JTextField();
        thresholdButton = new javax.swing.JButton();
        tpEvaluePanel = new javax.swing.JPanel();
        tpQvaluePanel = new javax.swing.JPanel();
        jSeparator4 = new javax.swing.JSeparator();
        siiComboBox = new javax.swing.JComboBox();
        siiListLabel = new javax.swing.JLabel();
        protocolPanel = new javax.swing.JPanel();
        protocolSummaryPanel = new javax.swing.JPanel();
        jMenuBar = new javax.swing.JMenuBar();
        fileMenu = new javax.swing.JMenu();
        openMenuItem = new javax.swing.JMenuItem();
        saveMenuItem = new javax.swing.JMenuItem();
        jSeparator3 = new javax.swing.JPopupMenu.Separator();
        exitMenuItem = new javax.swing.JMenuItem();
        exportMenu = new javax.swing.JMenu();
        exportProteinAmbiguityGrouMenuItem = new javax.swing.JMenuItem();
        exportProteinDetectionHypothesisMenuItem = new javax.swing.JMenuItem();
        exportPeptideSpectrumMatchesMenuItem = new javax.swing.JMenuItem();
        exportSeparator1 = new javax.swing.JPopupMenu.Separator();
        exportSpectrumIdentificationResultMenuItem = new javax.swing.JMenuItem();
        exportFragmentationMenuItem = new javax.swing.JMenuItem();
        exportSeparator2 = new javax.swing.JPopupMenu.Separator();
        exportFDR = new javax.swing.JMenuItem();
        optionsMenu = new javax.swing.JMenu();
        lookAndFeelMenu = new javax.swing.JMenu();
        lookAndFeelDefaultMenuItem = new javax.swing.JRadioButtonMenuItem();
        lookAndFeelMetalMenuItem = new javax.swing.JRadioButtonMenuItem();
        lookAndFeeMotifMenuItem = new javax.swing.JRadioButtonMenuItem();
        fontMenu = new javax.swing.JMenu();
        font12 = new javax.swing.JRadioButtonMenuItem();
        font14 = new javax.swing.JRadioButtonMenuItem();
        font16 = new javax.swing.JMenuItem();
        helpMenu = new javax.swing.JMenu();
        aboutMenuItem = new javax.swing.JMenuItem();

        fileChooser.setCurrentDirectory(null);
        fileChooser.addChoosableFileFilter (mzIdentMLFilter);

        setDefaultCloseOperation(javax.swing.WindowConstants.EXIT_ON_CLOSE);
        setTitle("MzIdentML Viewer");
        setMinimumSize(new java.awt.Dimension(900, 800));

        mainTabbedPane.setTabLayoutPolicy(javax.swing.JTabbedPane.SCROLL_TAB_LAYOUT);
        mainTabbedPane.setToolTipText("Global Statistics");
        mainTabbedPane.setPreferredSize(new java.awt.Dimension(889, 939));
        mainTabbedPane.addMouseListener(new java.awt.event.MouseAdapter() {
            public void mouseClicked(java.awt.event.MouseEvent evt) {
                mainTabbedPaneMouseClicked(evt);
            }
        });

        proteinViewPanel.setToolTipText("Protein View");
        proteinViewPanel.setName("Protein View"); // NOI18N
        proteinViewPanel.setPreferredSize(new java.awt.Dimension(889, 939));

        jProteinAmbiguityGroupPanel.setBorder(javax.swing.BorderFactory.createTitledBorder("Protein Group"));
        jProteinAmbiguityGroupPanel.setToolTipText("groups of proteins sharing some or all of the same peptides");
        jProteinAmbiguityGroupPanel.setPreferredSize(new java.awt.Dimension(772, 150));

        javax.swing.GroupLayout jProteinAmbiguityGroupPanelLayout = new javax.swing.GroupLayout(jProteinAmbiguityGroupPanel);
        jProteinAmbiguityGroupPanel.setLayout(jProteinAmbiguityGroupPanelLayout);
        jProteinAmbiguityGroupPanelLayout.setHorizontalGroup(
            jProteinAmbiguityGroupPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 358, Short.MAX_VALUE)
        );
        jProteinAmbiguityGroupPanelLayout.setVerticalGroup(
            jProteinAmbiguityGroupPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 213, Short.MAX_VALUE)
        );

        jProteinDetectionHypothesisPanel.setBorder(javax.swing.BorderFactory.createTitledBorder("Protein"));
        jProteinDetectionHypothesisPanel.setToolTipText("proteins inferred based on a set of peptide spectrum matches");
        jProteinDetectionHypothesisPanel.setPreferredSize(new java.awt.Dimension(772, 150));

        javax.swing.GroupLayout jProteinDetectionHypothesisPanelLayout = new javax.swing.GroupLayout(jProteinDetectionHypothesisPanel);
        jProteinDetectionHypothesisPanel.setLayout(jProteinDetectionHypothesisPanelLayout);
        jProteinDetectionHypothesisPanelLayout.setHorizontalGroup(
            jProteinDetectionHypothesisPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 351, Short.MAX_VALUE)
        );
        jProteinDetectionHypothesisPanelLayout.setVerticalGroup(
            jProteinDetectionHypothesisPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 608, Short.MAX_VALUE)
        );

        jProteinInfoPanel.setBorder(javax.swing.BorderFactory.createTitledBorder("Protein Info"));
        jProteinInfoPanel.setInheritsPopupMenu(true);

        jProteinSequencePanel.setBorder(javax.swing.BorderFactory.createTitledBorder("Protein Sequence"));

        jProteinSequenceTextPane.setContentType("text/html");
        jProteinSequenceScrollPane.setViewportView(jProteinSequenceTextPane);

        javax.swing.GroupLayout jProteinSequencePanelLayout = new javax.swing.GroupLayout(jProteinSequencePanel);
        jProteinSequencePanel.setLayout(jProteinSequencePanelLayout);
        jProteinSequencePanelLayout.setHorizontalGroup(
            jProteinSequencePanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addComponent(jProteinSequenceScrollPane, javax.swing.GroupLayout.DEFAULT_SIZE, 464, Short.MAX_VALUE)
        );
        jProteinSequencePanelLayout.setVerticalGroup(
            jProteinSequencePanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addComponent(jProteinSequenceScrollPane, javax.swing.GroupLayout.DEFAULT_SIZE, 159, Short.MAX_VALUE)
        );

        jSpectrumIdentificationItemProteinPanel.setBorder(javax.swing.BorderFactory.createTitledBorder("Peptide-Spectrum matches"));
        jSpectrumIdentificationItemProteinPanel.setToolTipText("Protein Detection Hypothesis");
        jSpectrumIdentificationItemProteinPanel.setPreferredSize(new java.awt.Dimension(772, 150));

        javax.swing.GroupLayout jSpectrumIdentificationItemProteinPanelLayout = new javax.swing.GroupLayout(jSpectrumIdentificationItemProteinPanel);
        jSpectrumIdentificationItemProteinPanel.setLayout(jSpectrumIdentificationItemProteinPanelLayout);
        jSpectrumIdentificationItemProteinPanelLayout.setHorizontalGroup(
            jSpectrumIdentificationItemProteinPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 464, Short.MAX_VALUE)
        );
        jSpectrumIdentificationItemProteinPanelLayout.setVerticalGroup(
            jSpectrumIdentificationItemProteinPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 514, Short.MAX_VALUE)
        );

        jScientificNameLabel.setText("Scientific name:");

        jProteinDescriptionPanel.setBorder(javax.swing.BorderFactory.createTitledBorder("Protein Description"));

        jProteinDescriptionEditorPane.setContentType("text/html");
        jProteinDescriptionScrollPane.setViewportView(jProteinDescriptionEditorPane);

        javax.swing.GroupLayout jProteinDescriptionPanelLayout = new javax.swing.GroupLayout(jProteinDescriptionPanel);
        jProteinDescriptionPanel.setLayout(jProteinDescriptionPanelLayout);
        jProteinDescriptionPanelLayout.setHorizontalGroup(
            jProteinDescriptionPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addComponent(jProteinDescriptionScrollPane, javax.swing.GroupLayout.DEFAULT_SIZE, 464, Short.MAX_VALUE)
        );
        jProteinDescriptionPanelLayout.setVerticalGroup(
            jProteinDescriptionPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addComponent(jProteinDescriptionScrollPane, javax.swing.GroupLayout.DEFAULT_SIZE, 59, Short.MAX_VALUE)
        );

        javax.swing.GroupLayout jProteinInfoPanelLayout = new javax.swing.GroupLayout(jProteinInfoPanel);
        jProteinInfoPanel.setLayout(jProteinInfoPanelLayout);
        jProteinInfoPanelLayout.setHorizontalGroup(
            jProteinInfoPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jProteinInfoPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addComponent(jScientificNameLabel)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addComponent(jScientificNameValueLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 157, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addContainerGap(231, Short.MAX_VALUE))
            .addComponent(jProteinDescriptionPanel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
            .addComponent(jProteinSequencePanel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
            .addComponent(jSpectrumIdentificationItemProteinPanel, javax.swing.GroupLayout.DEFAULT_SIZE, 483, Short.MAX_VALUE)
        );
        jProteinInfoPanelLayout.setVerticalGroup(
            jProteinInfoPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jProteinInfoPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(jProteinInfoPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jScientificNameLabel)
                    .addComponent(jScientificNameValueLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 14, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addGap(11, 11, 11)
                .addComponent(jProteinDescriptionPanel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addComponent(jProteinSequencePanel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addGap(18, 18, 18)
                .addComponent(jSpectrumIdentificationItemProteinPanel, javax.swing.GroupLayout.DEFAULT_SIZE, 537, Short.MAX_VALUE))
        );

        jSpectrumIdentificationItemProteinPanel.getAccessibleContext().setAccessibleDescription("Spectrum Identification Item");

        javax.swing.GroupLayout proteinViewPanelLayout = new javax.swing.GroupLayout(proteinViewPanel);
        proteinViewPanel.setLayout(proteinViewPanelLayout);
        proteinViewPanelLayout.setHorizontalGroup(
            proteinViewPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(proteinViewPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(proteinViewPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(jProteinDetectionHypothesisPanel, javax.swing.GroupLayout.DEFAULT_SIZE, 363, Short.MAX_VALUE)
                    .addComponent(jProteinAmbiguityGroupPanel, javax.swing.GroupLayout.DEFAULT_SIZE, 363, Short.MAX_VALUE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(jProteinInfoPanel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addContainerGap())
        );
        proteinViewPanelLayout.setVerticalGroup(
            proteinViewPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, proteinViewPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(proteinViewPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                    .addComponent(jProteinInfoPanel, javax.swing.GroupLayout.Alignment.LEADING, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addGroup(proteinViewPanelLayout.createSequentialGroup()
                        .addComponent(jProteinAmbiguityGroupPanel, javax.swing.GroupLayout.PREFERRED_SIZE, 236, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addGap(18, 18, 18)
                        .addComponent(jProteinDetectionHypothesisPanel, javax.swing.GroupLayout.DEFAULT_SIZE, 635, Short.MAX_VALUE)))
                .addContainerGap())
        );

        mainTabbedPane.addTab("Protein View", null, proteinViewPanel, "Protein View");
        proteinViewPanel.getAccessibleContext().setAccessibleName("Protein View");
        proteinViewPanel.getAccessibleContext().setAccessibleParent(mainTabbedPane);

        spectrumViewPanel.setToolTipText("Spectrum View");
        spectrumViewPanel.setPreferredSize(new java.awt.Dimension(889, 939));

        jSpectrumIdentificationResultPanel.setBorder(javax.swing.BorderFactory.createTitledBorder("Spectrum List"));
        jSpectrumIdentificationResultPanel.setToolTipText("Protein Ambiguity Group");
        jSpectrumIdentificationResultPanel.setMinimumSize(new java.awt.Dimension(404, 569));

        javax.swing.GroupLayout jSpectrumIdentificationResultPanelLayout = new javax.swing.GroupLayout(jSpectrumIdentificationResultPanel);
        jSpectrumIdentificationResultPanel.setLayout(jSpectrumIdentificationResultPanelLayout);
        jSpectrumIdentificationResultPanelLayout.setHorizontalGroup(
            jSpectrumIdentificationResultPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 443, Short.MAX_VALUE)
        );
        jSpectrumIdentificationResultPanelLayout.setVerticalGroup(
            jSpectrumIdentificationResultPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 0, Short.MAX_VALUE)
        );

        jSpectrumIdentificationItemPanel.setBorder(javax.swing.BorderFactory.createTitledBorder("Peptide-Spectrum matches"));
        jSpectrumIdentificationItemPanel.setToolTipText("Spectrum Identification Item");
        jSpectrumIdentificationItemPanel.setAutoscrolls(true);

        javax.swing.GroupLayout jSpectrumIdentificationItemPanelLayout = new javax.swing.GroupLayout(jSpectrumIdentificationItemPanel);
        jSpectrumIdentificationItemPanel.setLayout(jSpectrumIdentificationItemPanelLayout);
        jSpectrumIdentificationItemPanelLayout.setHorizontalGroup(
            jSpectrumIdentificationItemPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 444, Short.MAX_VALUE)
        );
        jSpectrumIdentificationItemPanelLayout.setVerticalGroup(
            jSpectrumIdentificationItemPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 418, Short.MAX_VALUE)
        );

        jPeptideEvidencePanel.setBorder(javax.swing.BorderFactory.createTitledBorder("Peptide Evidence"));
        jPeptideEvidencePanel.setToolTipText("Peptide Evidence");
        jPeptideEvidencePanel.setAutoscrolls(true);

        javax.swing.GroupLayout jPeptideEvidencePanelLayout = new javax.swing.GroupLayout(jPeptideEvidencePanel);
        jPeptideEvidencePanel.setLayout(jPeptideEvidencePanelLayout);
        jPeptideEvidencePanelLayout.setHorizontalGroup(
            jPeptideEvidencePanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 443, Short.MAX_VALUE)
        );
        jPeptideEvidencePanelLayout.setVerticalGroup(
            jPeptideEvidencePanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 177, Short.MAX_VALUE)
        );

        jSpectrumPanel.setBorder(javax.swing.BorderFactory.createTitledBorder("Spectrum"));
        jSpectrumPanel.setAutoscrolls(true);
        jSpectrumPanel.setPreferredSize(new java.awt.Dimension(362, 569));

        jFragmentationPanel.setBorder(javax.swing.BorderFactory.createTitledBorder("Fragmentation"));
        jFragmentationPanel.setAutoscrolls(true);
        jFragmentationPanel.setMinimumSize(new java.awt.Dimension(0, 0));
        jFragmentationPanel.setPreferredSize(new java.awt.Dimension(383, 447));

        javax.swing.GroupLayout jFragmentationPanelLayout = new javax.swing.GroupLayout(jFragmentationPanel);
        jFragmentationPanel.setLayout(jFragmentationPanelLayout);
        jFragmentationPanelLayout.setHorizontalGroup(
            jFragmentationPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 379, Short.MAX_VALUE)
        );
        jFragmentationPanelLayout.setVerticalGroup(
            jFragmentationPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 218, Short.MAX_VALUE)
        );

        jGraph.setBorder(javax.swing.BorderFactory.createTitledBorder("Graph"));

        javax.swing.GroupLayout jGraphLayout = new javax.swing.GroupLayout(jGraph);
        jGraph.setLayout(jGraphLayout);
        jGraphLayout.setHorizontalGroup(
            jGraphLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 379, Short.MAX_VALUE)
        );
        jGraphLayout.setVerticalGroup(
            jGraphLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 411, Short.MAX_VALUE)
        );

        jExperimentalFilterPanel.setBorder(javax.swing.BorderFactory.createTitledBorder("Experimental Filtering"));

        javax.swing.GroupLayout jExperimentalFilterPanelLayout = new javax.swing.GroupLayout(jExperimentalFilterPanel);
        jExperimentalFilterPanel.setLayout(jExperimentalFilterPanelLayout);
        jExperimentalFilterPanelLayout.setHorizontalGroup(
            jExperimentalFilterPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 379, Short.MAX_VALUE)
        );
        jExperimentalFilterPanelLayout.setVerticalGroup(
            jExperimentalFilterPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 30, Short.MAX_VALUE)
        );

        jTheoreticalFilterPanel.setBorder(javax.swing.BorderFactory.createTitledBorder("Theoretical Filtering"));

        bCheckBox.setText("b");
        bCheckBox.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                bCheckBoxActionPerformed(evt);
            }
        });

        yCheckBox.setText("y");
        yCheckBox.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                yCheckBoxActionPerformed(evt);
            }
        });

        oneCheckBox.setText("1");
        oneCheckBox.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                oneCheckBoxActionPerformed(evt);
            }
        });

        twoCheckBox.setText("2");
        twoCheckBox.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                twoCheckBoxActionPerformed(evt);
            }
        });

        javax.swing.GroupLayout jTheoreticalFilterPanelLayout = new javax.swing.GroupLayout(jTheoreticalFilterPanel);
        jTheoreticalFilterPanel.setLayout(jTheoreticalFilterPanelLayout);
        jTheoreticalFilterPanelLayout.setHorizontalGroup(
            jTheoreticalFilterPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jTheoreticalFilterPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addComponent(bCheckBox)
                .addGap(64, 64, 64)
                .addComponent(yCheckBox)
                .addGap(61, 61, 61)
                .addComponent(oneCheckBox)
                .addGap(56, 56, 56)
                .addComponent(twoCheckBox)
                .addContainerGap(67, Short.MAX_VALUE))
        );
        jTheoreticalFilterPanelLayout.setVerticalGroup(
            jTheoreticalFilterPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jTheoreticalFilterPanelLayout.createSequentialGroup()
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addGroup(jTheoreticalFilterPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(bCheckBox)
                    .addComponent(yCheckBox)
                    .addComponent(oneCheckBox)
                    .addComponent(twoCheckBox)))
        );

        javax.swing.GroupLayout jSpectrumPanelLayout = new javax.swing.GroupLayout(jSpectrumPanel);
        jSpectrumPanel.setLayout(jSpectrumPanelLayout);
        jSpectrumPanelLayout.setHorizontalGroup(
            jSpectrumPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addComponent(jGraph, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
            .addComponent(jExperimentalFilterPanel, javax.swing.GroupLayout.Alignment.TRAILING, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
            .addComponent(jTheoreticalFilterPanel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
            .addComponent(jFragmentationPanel, javax.swing.GroupLayout.DEFAULT_SIZE, 390, Short.MAX_VALUE)
        );
        jSpectrumPanelLayout.setVerticalGroup(
            jSpectrumPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jSpectrumPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addComponent(jGraph, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addGap(18, 18, 18)
                .addComponent(jExperimentalFilterPanel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addGap(18, 18, 18)
                .addComponent(jTheoreticalFilterPanel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addGap(18, 18, 18)
                .addComponent(jFragmentationPanel, javax.swing.GroupLayout.DEFAULT_SIZE, 261, Short.MAX_VALUE))
        );

        javax.swing.GroupLayout spectrumViewPanelLayout = new javax.swing.GroupLayout(spectrumViewPanel);
        spectrumViewPanel.setLayout(spectrumViewPanelLayout);
        spectrumViewPanelLayout.setHorizontalGroup(
            spectrumViewPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(spectrumViewPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(spectrumViewPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(jSpectrumIdentificationItemPanel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(jSpectrumIdentificationResultPanel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(jPeptideEvidencePanel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(jSpectrumPanel, javax.swing.GroupLayout.DEFAULT_SIZE, 402, Short.MAX_VALUE)
                .addContainerGap())
        );
        spectrumViewPanelLayout.setVerticalGroup(
            spectrumViewPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(spectrumViewPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(spectrumViewPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(spectrumViewPanelLayout.createSequentialGroup()
                        .addComponent(jSpectrumIdentificationResultPanel, javax.swing.GroupLayout.PREFERRED_SIZE, 226, Short.MAX_VALUE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                        .addComponent(jSpectrumIdentificationItemPanel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                        .addComponent(jPeptideEvidencePanel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                    .addComponent(jSpectrumPanel, javax.swing.GroupLayout.Alignment.TRAILING, javax.swing.GroupLayout.DEFAULT_SIZE, 889, Short.MAX_VALUE))
                .addContainerGap())
        );

        jSpectrumIdentificationResultPanel.getAccessibleContext().setAccessibleDescription("Spectrum Identification Result");

        mainTabbedPane.addTab("Spectrum Summary", null, spectrumViewPanel, "Spectrum Summary");
        spectrumViewPanel.getAccessibleContext().setAccessibleName("Spectrum Summary");
        spectrumViewPanel.getAccessibleContext().setAccessibleDescription("Spectrum Summary");
        spectrumViewPanel.getAccessibleContext().setAccessibleParent(mainTabbedPane);

        proteinDBViewPanel.setPreferredSize(new java.awt.Dimension(889, 939));

        dBSequencePanel.setBorder(javax.swing.BorderFactory.createTitledBorder("DB Sequence"));

        javax.swing.GroupLayout dBSequencePanelLayout = new javax.swing.GroupLayout(dBSequencePanel);
        dBSequencePanel.setLayout(dBSequencePanelLayout);
        dBSequencePanelLayout.setHorizontalGroup(
            dBSequencePanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 852, Short.MAX_VALUE)
        );
        dBSequencePanelLayout.setVerticalGroup(
            dBSequencePanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 866, Short.MAX_VALUE)
        );

        javax.swing.GroupLayout proteinDBViewPanelLayout = new javax.swing.GroupLayout(proteinDBViewPanel);
        proteinDBViewPanel.setLayout(proteinDBViewPanelLayout);
        proteinDBViewPanelLayout.setHorizontalGroup(
            proteinDBViewPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(proteinDBViewPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addComponent(dBSequencePanel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addContainerGap())
        );
        proteinDBViewPanelLayout.setVerticalGroup(
            proteinDBViewPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(proteinDBViewPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addComponent(dBSequencePanel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addContainerGap())
        );

        mainTabbedPane.addTab("Protein DB View", proteinDBViewPanel);
        proteinDBViewPanel.getAccessibleContext().setAccessibleParent(mainTabbedPane);

        globalStatisticsPanel.setPreferredSize(new java.awt.Dimension(889, 939));

        summaryPanel.setBorder(javax.swing.BorderFactory.createTitledBorder("Summary"));

        totalSIRLabel.setText("Total SIR:");

        totalSIRLabelValue.setText("0");

        totalSIILabel.setText("Total PSM:");

        totalSIIbelowThresholdLabel.setText("Total PSM below Threshold:");

        totalSIIbelowThresholdLabelValue.setText("0");

        totalSIIaboveThresholdLabel.setText("Total PSM pass Threshold:");

        totalSIIaboveThresholdLabelValue.setText("0");

        totalSIIaboveThresholdRankOneLabel.setText("Total PSM pass Threshold and where rank=1:");

        totalSIIaboveThresholdRankOneLabelValue.setText("0");

        percentIdentifiedSpectraLabel.setText("Percent identified spectra:");

        percentIdentifiedSpectraLabelValue.setText("0");

        totalPeptidesaboveThresholdLabel.setText("Total non-redundant peptides above Threshold:");

        totalPeptidesaboveThresholdLabelValue.setText("0");

        totalPAGsLabel.setText("Total PAGs:");

        totalPAGsLabelValue.setText("0");

        totalPDHsLabel.setText("Total PDHs:");

        totalPDHsaboveThresholdLabel.setText("Total PDH above Threshold:");

        totalPDHsaboveThresholdLabelValue.setText("0");

        isDecoySii.setText("PSM with Decoy = true:");

        isDecoySiiValue.setText("0");

        isDecoySiiFalse.setText("PSM with Decoy = false:");

        isDecoySiiFalseValue.setText("0");

        manualDecoyLabel.setText("Manual Decoy:");

        manualDecoyPrefix.setText("Prefix");

        manualDecoyRatio.setText("Ratio");

        manualDecoyPrefixValue.setText("Rev_");
        manualDecoyPrefixValue.setEnabled(false);

        manualDecoyRatioValue.setText("1");
        manualDecoyRatioValue.setEnabled(false);

        manualDecoy.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                manualDecoyActionPerformed(evt);
            }
        });

        fpSiiLabel.setText("FP for PSM:");

        fpSiiValue.setText("0");

        tpSiiLabel.setText("TP for PSM:");

        tpSiiValue.setText("0");

        fdrSiiLabel.setText("FDR for PSM:");

        fdrSiiValue.setText("0");

        fpProteinsLabel.setText("FP for proteins:");

        tpProteinsLabel.setText("TP for proteins:");

        fdrProteinsLabel.setText("FDR for proteins:");

        fpProteinsValue.setText("0");

        tpProteinsValue.setText("0");

        fdrProteinsValue.setText("0");

        manualCalculate.setText("Calculate");
        manualCalculate.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                manualCalculateActionPerformed(evt);
            }
        });

        fdrPanel.setBorder(javax.swing.BorderFactory.createTitledBorder("FDR Graph"));
        fdrPanel.setPreferredSize(new java.awt.Dimension(257, 255));

        javax.swing.GroupLayout fdrPanelLayout = new javax.swing.GroupLayout(fdrPanel);
        fdrPanel.setLayout(fdrPanelLayout);
        fdrPanelLayout.setHorizontalGroup(
            fdrPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 265, Short.MAX_VALUE)
        );
        fdrPanelLayout.setVerticalGroup(
            fdrPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 389, Short.MAX_VALUE)
        );

        thresholdLabel.setText("Threshold");

        thresholdValue.setEnabled(false);

        thresholdButton.setText("Change");
        thresholdButton.setEnabled(false);

        tpEvaluePanel.setBorder(javax.swing.BorderFactory.createTitledBorder("TP vs FP vs E-value"));
        tpEvaluePanel.setPreferredSize(new java.awt.Dimension(257, 255));

        javax.swing.GroupLayout tpEvaluePanelLayout = new javax.swing.GroupLayout(tpEvaluePanel);
        tpEvaluePanel.setLayout(tpEvaluePanelLayout);
        tpEvaluePanelLayout.setHorizontalGroup(
            tpEvaluePanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 261, Short.MAX_VALUE)
        );
        tpEvaluePanelLayout.setVerticalGroup(
            tpEvaluePanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 397, Short.MAX_VALUE)
        );

        tpQvaluePanel.setBorder(javax.swing.BorderFactory.createTitledBorder("TP vs Q-value"));
        tpQvaluePanel.setPreferredSize(new java.awt.Dimension(257, 255));

        javax.swing.GroupLayout tpQvaluePanelLayout = new javax.swing.GroupLayout(tpQvaluePanel);
        tpQvaluePanel.setLayout(tpQvaluePanelLayout);
        tpQvaluePanelLayout.setHorizontalGroup(
            tpQvaluePanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 250, Short.MAX_VALUE)
        );
        tpQvaluePanelLayout.setVerticalGroup(
            tpQvaluePanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 389, Short.MAX_VALUE)
        );

        jSeparator4.setOrientation(javax.swing.SwingConstants.VERTICAL);

        siiComboBox.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                siiComboBoxActionPerformed(evt);
            }
        });

        siiListLabel.setText("SII List:");

        javax.swing.GroupLayout summaryPanelLayout = new javax.swing.GroupLayout(summaryPanel);
        summaryPanel.setLayout(summaryPanelLayout);
        summaryPanelLayout.setHorizontalGroup(
            summaryPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addComponent(jSeparator2, javax.swing.GroupLayout.DEFAULT_SIZE, 852, Short.MAX_VALUE)
            .addComponent(jSeparator1, javax.swing.GroupLayout.Alignment.TRAILING, javax.swing.GroupLayout.DEFAULT_SIZE, 852, Short.MAX_VALUE)
            .addGroup(summaryPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(summaryPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(summaryPanelLayout.createSequentialGroup()
                        .addComponent(fdrPanel, javax.swing.GroupLayout.DEFAULT_SIZE, 275, Short.MAX_VALUE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                        .addComponent(tpEvaluePanel, javax.swing.GroupLayout.DEFAULT_SIZE, 273, Short.MAX_VALUE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                        .addComponent(tpQvaluePanel, javax.swing.GroupLayout.DEFAULT_SIZE, 264, Short.MAX_VALUE))
                    .addGroup(summaryPanelLayout.createSequentialGroup()
                        .addGroup(summaryPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addGroup(summaryPanelLayout.createSequentialGroup()
                                .addGroup(summaryPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                    .addComponent(isDecoySii)
                                    .addComponent(isDecoySiiFalse)
                                    .addComponent(fpSiiLabel)
                                    .addComponent(tpSiiLabel)
                                    .addComponent(fdrSiiLabel)
                                    .addGroup(summaryPanelLayout.createSequentialGroup()
                                        .addComponent(manualDecoy)
                                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                                        .addComponent(manualDecoyLabel))
                                    .addComponent(fpProteinsLabel)
                                    .addComponent(tpProteinsLabel)
                                    .addComponent(fdrProteinsLabel)
                                    .addComponent(thresholdLabel))
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addGroup(summaryPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                    .addComponent(fdrProteinsValue)
                                    .addComponent(tpProteinsValue)
                                    .addComponent(fpProteinsValue)
                                    .addGroup(summaryPanelLayout.createSequentialGroup()
                                        .addGroup(summaryPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                            .addComponent(isDecoySiiValue)
                                            .addComponent(isDecoySiiFalseValue)
                                            .addComponent(manualDecoyPrefix)
                                            .addComponent(fpSiiValue)
                                            .addComponent(fdrSiiValue)
                                            .addGroup(summaryPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING, false)
                                                .addComponent(tpSiiValue, javax.swing.GroupLayout.Alignment.LEADING, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                                .addComponent(manualDecoyRatio, javax.swing.GroupLayout.Alignment.LEADING, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)))
                                        .addGap(212, 212, 212)
                                        .addGroup(summaryPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                                            .addComponent(manualDecoyRatioValue, javax.swing.GroupLayout.DEFAULT_SIZE, 51, Short.MAX_VALUE)
                                            .addComponent(manualDecoyPrefixValue))
                                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                        .addComponent(manualCalculate))
                                    .addGroup(summaryPanelLayout.createSequentialGroup()
                                        .addComponent(thresholdValue, javax.swing.GroupLayout.PREFERRED_SIZE, 53, javax.swing.GroupLayout.PREFERRED_SIZE)
                                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                        .addComponent(thresholdButton))))
                            .addGroup(summaryPanelLayout.createSequentialGroup()
                                .addGroup(summaryPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                    .addComponent(totalSIRLabel)
                                    .addComponent(totalSIILabel)
                                    .addComponent(siiListLabel))
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addGroup(summaryPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                    .addGroup(summaryPanelLayout.createSequentialGroup()
                                        .addGroup(summaryPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                            .addComponent(totalSIRLabelValue)
                                            .addComponent(totalSIIbelowThresholdLabel)
                                            .addComponent(totalSIIaboveThresholdLabel)
                                            .addComponent(totalSIIaboveThresholdRankOneLabel)
                                            .addComponent(percentIdentifiedSpectraLabel)
                                            .addComponent(totalPeptidesaboveThresholdLabel))
                                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                                        .addGroup(summaryPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                            .addComponent(totalPeptidesaboveThresholdLabelValue)
                                            .addComponent(percentIdentifiedSpectraLabelValue)
                                            .addComponent(totalSIIaboveThresholdRankOneLabelValue)
                                            .addComponent(totalSIIaboveThresholdLabelValue)
                                            .addComponent(totalSIIbelowThresholdLabelValue)))
                                    .addComponent(siiComboBox, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                                .addGap(25, 25, 25)
                                .addComponent(jSeparator4, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                                .addGap(28, 28, 28)
                                .addGroup(summaryPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                    .addComponent(totalPAGsLabel)
                                    .addComponent(totalPDHsLabel))
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addGroup(summaryPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                    .addComponent(totalPAGsLabelValue)
                                    .addGroup(summaryPanelLayout.createSequentialGroup()
                                        .addComponent(totalPDHsaboveThresholdLabel)
                                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                                        .addComponent(totalPDHsaboveThresholdLabelValue)))))
                        .addGap(0, 258, Short.MAX_VALUE)))
                .addContainerGap())
        );
        summaryPanelLayout.setVerticalGroup(
            summaryPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(summaryPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(summaryPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(summaryPanelLayout.createSequentialGroup()
                        .addGroup(summaryPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                            .addComponent(siiComboBox, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addComponent(siiListLabel))
                        .addGap(8, 8, 8)
                        .addGroup(summaryPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                            .addGroup(summaryPanelLayout.createSequentialGroup()
                                .addGroup(summaryPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                                    .addComponent(totalSIRLabel)
                                    .addComponent(totalSIRLabelValue))
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addGroup(summaryPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                                    .addComponent(totalSIILabel)
                                    .addComponent(totalSIIbelowThresholdLabel)
                                    .addComponent(totalSIIbelowThresholdLabelValue)))
                            .addGroup(summaryPanelLayout.createSequentialGroup()
                                .addGroup(summaryPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                                    .addComponent(totalPAGsLabel)
                                    .addComponent(totalPAGsLabelValue))
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addGroup(summaryPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                                    .addComponent(totalPDHsLabel)
                                    .addComponent(totalPDHsaboveThresholdLabel)
                                    .addComponent(totalPDHsaboveThresholdLabelValue))))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(summaryPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                            .addComponent(totalSIIaboveThresholdLabel)
                            .addComponent(totalSIIaboveThresholdLabelValue))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(summaryPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                            .addComponent(totalSIIaboveThresholdRankOneLabel)
                            .addComponent(totalSIIaboveThresholdRankOneLabelValue))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(summaryPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                            .addComponent(percentIdentifiedSpectraLabelValue)
                            .addComponent(percentIdentifiedSpectraLabel))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(summaryPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                            .addComponent(totalPeptidesaboveThresholdLabel)
                            .addComponent(totalPeptidesaboveThresholdLabelValue)))
                    .addComponent(jSeparator4, javax.swing.GroupLayout.PREFERRED_SIZE, 142, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addGap(18, 18, 18)
                .addComponent(jSeparator1, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(summaryPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(thresholdLabel)
                    .addComponent(thresholdValue, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(thresholdButton))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(summaryPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(isDecoySii)
                    .addComponent(isDecoySiiValue))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(summaryPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(isDecoySiiFalse)
                    .addComponent(isDecoySiiFalseValue))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(summaryPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                    .addGroup(summaryPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                        .addComponent(manualDecoyPrefix)
                        .addComponent(manualDecoyPrefixValue, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addComponent(manualDecoyLabel))
                    .addComponent(manualDecoy))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(summaryPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(manualDecoyRatio)
                    .addComponent(manualDecoyRatioValue, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(manualCalculate))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(summaryPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(fpSiiLabel)
                    .addComponent(fpSiiValue))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(summaryPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(tpSiiLabel)
                    .addComponent(tpSiiValue))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(summaryPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(fdrSiiLabel)
                    .addComponent(fdrSiiValue))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(summaryPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(fpProteinsLabel)
                    .addComponent(fpProteinsValue))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(summaryPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(tpProteinsLabel)
                    .addComponent(tpProteinsValue))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(summaryPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(fdrProteinsLabel)
                    .addComponent(fdrProteinsValue))
                .addGap(10, 10, 10)
                .addComponent(jSeparator2, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(summaryPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                    .addComponent(tpEvaluePanel, javax.swing.GroupLayout.DEFAULT_SIZE, 420, Short.MAX_VALUE)
                    .addComponent(fdrPanel, javax.swing.GroupLayout.DEFAULT_SIZE, 420, Short.MAX_VALUE)
                    .addComponent(tpQvaluePanel, javax.swing.GroupLayout.DEFAULT_SIZE, 420, Short.MAX_VALUE))
                .addContainerGap())
        );

        javax.swing.GroupLayout globalStatisticsPanelLayout = new javax.swing.GroupLayout(globalStatisticsPanel);
        globalStatisticsPanel.setLayout(globalStatisticsPanelLayout);
        globalStatisticsPanelLayout.setHorizontalGroup(
            globalStatisticsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(globalStatisticsPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addComponent(summaryPanel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addContainerGap())
        );
        globalStatisticsPanelLayout.setVerticalGroup(
            globalStatisticsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(globalStatisticsPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addComponent(summaryPanel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addContainerGap())
        );

        mainTabbedPane.addTab("Global Statistics", null, globalStatisticsPanel, "Global Statistics");
        globalStatisticsPanel.getAccessibleContext().setAccessibleName("Global Statistics");
        globalStatisticsPanel.getAccessibleContext().setAccessibleDescription("Global Statistics");
        globalStatisticsPanel.getAccessibleContext().setAccessibleParent(mainTabbedPane);

        protocolPanel.setPreferredSize(new java.awt.Dimension(889, 939));

        protocolSummaryPanel.setBorder(javax.swing.BorderFactory.createTitledBorder("Summary"));

        javax.swing.GroupLayout protocolSummaryPanelLayout = new javax.swing.GroupLayout(protocolSummaryPanel);
        protocolSummaryPanel.setLayout(protocolSummaryPanelLayout);
        protocolSummaryPanelLayout.setHorizontalGroup(
            protocolSummaryPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 852, Short.MAX_VALUE)
        );
        protocolSummaryPanelLayout.setVerticalGroup(
            protocolSummaryPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 866, Short.MAX_VALUE)
        );

        javax.swing.GroupLayout protocolPanelLayout = new javax.swing.GroupLayout(protocolPanel);
        protocolPanel.setLayout(protocolPanelLayout);
        protocolPanelLayout.setHorizontalGroup(
            protocolPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(protocolPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addComponent(protocolSummaryPanel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addContainerGap())
        );
        protocolPanelLayout.setVerticalGroup(
            protocolPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(protocolPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addComponent(protocolSummaryPanel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addContainerGap())
        );

        mainTabbedPane.addTab("Protocols", protocolPanel);
        protocolPanel.getAccessibleContext().setAccessibleParent(mainTabbedPane);

        javax.swing.GroupLayout mainPanelLayout = new javax.swing.GroupLayout(mainPanel);
        mainPanel.setLayout(mainPanelLayout);
        mainPanelLayout.setHorizontalGroup(
            mainPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(mainPanelLayout.createSequentialGroup()
                .addComponent(mainTabbedPane, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addContainerGap())
        );
        mainPanelLayout.setVerticalGroup(
            mainPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(mainPanelLayout.createSequentialGroup()
                .addComponent(mainTabbedPane, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addContainerGap())
        );

        fileMenu.setText("File");

        openMenuItem.setAccelerator(javax.swing.KeyStroke.getKeyStroke(java.awt.event.KeyEvent.VK_O, java.awt.event.InputEvent.SHIFT_MASK | java.awt.event.InputEvent.CTRL_MASK));
        openMenuItem.setText("Open");
        openMenuItem.setToolTipText("Open");
        openMenuItem.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                openMenuItemActionPerformed(evt);
            }
        });
        fileMenu.add(openMenuItem);

        saveMenuItem.setAccelerator(javax.swing.KeyStroke.getKeyStroke(java.awt.event.KeyEvent.VK_S, java.awt.event.InputEvent.SHIFT_MASK | java.awt.event.InputEvent.CTRL_MASK));
        saveMenuItem.setText("Save");
        saveMenuItem.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                saveMenuItemActionPerformed(evt);
            }
        });
        fileMenu.add(saveMenuItem);
        fileMenu.add(jSeparator3);

        exitMenuItem.setAccelerator(javax.swing.KeyStroke.getKeyStroke(java.awt.event.KeyEvent.VK_Q, java.awt.event.InputEvent.CTRL_MASK));
        exitMenuItem.setText("Exit");
        exitMenuItem.setToolTipText("Exit");
        exitMenuItem.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                exitMenuItemActionPerformed(evt);
            }
        });
        fileMenu.add(exitMenuItem);

        jMenuBar.add(fileMenu);

        exportMenu.setText("Export");

        exportProteinAmbiguityGrouMenuItem.setText("Export Protein Ambiguity Groups");
        exportProteinAmbiguityGrouMenuItem.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                exportProteinAmbiguityGrouMenuItemActionPerformed(evt);
            }
        });
        exportMenu.add(exportProteinAmbiguityGrouMenuItem);

        exportProteinDetectionHypothesisMenuItem.setText("Export Protein Detection Hypothesis");
        exportProteinDetectionHypothesisMenuItem.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                exportProteinDetectionHypothesisMenuItemActionPerformed(evt);
            }
        });
        exportMenu.add(exportProteinDetectionHypothesisMenuItem);

        exportPeptideSpectrumMatchesMenuItem.setText("Export Peptide-Spectrum Matches");
        exportPeptideSpectrumMatchesMenuItem.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                exportPeptideSpectrumMatchesMenuItemActionPerformed(evt);
            }
        });
        exportMenu.add(exportPeptideSpectrumMatchesMenuItem);
        exportMenu.add(exportSeparator1);

        exportSpectrumIdentificationResultMenuItem.setText("Export Spectrum Identification Result");
        exportSpectrumIdentificationResultMenuItem.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                exportSpectrumIdentificationResultMenuItemActionPerformed(evt);
            }
        });
        exportMenu.add(exportSpectrumIdentificationResultMenuItem);

        exportFragmentationMenuItem.setText("Export Fragmentation");
        exportFragmentationMenuItem.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                exportFragmentationMenuItemActionPerformed(evt);
            }
        });
        exportMenu.add(exportFragmentationMenuItem);
        exportMenu.add(exportSeparator2);

        exportFDR.setText("Export FDR as CSV");
        exportFDR.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                exportFDRActionPerformed(evt);
            }
        });
        exportMenu.add(exportFDR);

        jMenuBar.add(exportMenu);

        optionsMenu.setText("Options");

        lookAndFeelMenu.setText("Look & Feel");

        lookAndFeelDefaultMenuItem.setSelected(true);
        lookAndFeelDefaultMenuItem.setText("Default");
        lookAndFeelDefaultMenuItem.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                lookAndFeelDefaultMenuItemActionPerformed(evt);
            }
        });
        lookAndFeelMenu.add(lookAndFeelDefaultMenuItem);

        lookAndFeelMetalMenuItem.setText("Metal");
        lookAndFeelMetalMenuItem.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                lookAndFeelMetalMenuItemActionPerformed(evt);
            }
        });
        lookAndFeelMenu.add(lookAndFeelMetalMenuItem);

        lookAndFeeMotifMenuItem.setText("Motif");
        lookAndFeeMotifMenuItem.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                lookAndFeeMotifMenuItemActionPerformed(evt);
            }
        });
        lookAndFeelMenu.add(lookAndFeeMotifMenuItem);

        optionsMenu.add(lookAndFeelMenu);

        fontMenu.setText("Font");

        font12.setText("12");
        font12.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                font12ActionPerformed(evt);
            }
        });
        fontMenu.add(font12);

        font14.setText("14");
        font14.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                font14ActionPerformed(evt);
            }
        });
        fontMenu.add(font14);

        font16.setText("16");
        font16.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                font16ActionPerformed(evt);
            }
        });
        fontMenu.add(font16);

        optionsMenu.add(fontMenu);

        jMenuBar.add(optionsMenu);

        helpMenu.setText("Help");

        aboutMenuItem.setAccelerator(javax.swing.KeyStroke.getKeyStroke(java.awt.event.KeyEvent.VK_A, java.awt.event.InputEvent.CTRL_MASK));
        aboutMenuItem.setText("About");
        aboutMenuItem.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                aboutMenuItemActionPerformed(evt);
            }
        });
        helpMenu.add(aboutMenuItem);

        jMenuBar.add(helpMenu);

        setJMenuBar(jMenuBar);

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(getContentPane());
        getContentPane().setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addComponent(mainPanel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addComponent(mainPanel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
        );

        getAccessibleContext().setAccessibleDescription("MzIdentML Viewer");

        pack();
    }// </editor-fold>//GEN-END:initComponents

    /**
     * exit Menu Item Action Performed
     */
    private void exitMenuItemActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_exitMenuItemActionPerformed
        int exit = JOptionPane.showConfirmDialog(this, "Are you sure?", "Close mzIdentML Viewer", JOptionPane.YES_NO_OPTION);
        if (exit == JOptionPane.YES_OPTION) {
            System.exit(0);
        }
    }//GEN-LAST:event_exitMenuItemActionPerformed
    /**
     * about Menu Item Action Performed
     */
    private void aboutMenuItemActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_aboutMenuItemActionPerformed
        JOptionPane.showMessageDialog(this, "mzIdentML Viewer", "About", JOptionPane.INFORMATION_MESSAGE);
    }//GEN-LAST:event_aboutMenuItemActionPerformed
    /**
     * protein Detection Hypothesis Table Mouse Clicked
     */
    private void proteinDetectionHypothesisTableMouseClicked(MouseEvent evt) {
        spectrumIdentificationItemProteinViewTable.setSortable(false);
        this.setCursor(new java.awt.Cursor(java.awt.Cursor.WAIT_CURSOR));
        int row = proteinDetectionHypothesisTable.getSelectedRow();
        SpectrumIdentificationItem spectrumIdentificationItem2 = null;
        if (row != -1) {
            try {
                while (spectrumIdentificationItemProteinViewTable.getRowCount() > 0) {
                    ((DefaultTableModel) spectrumIdentificationItemProteinViewTable.getModel()).removeRow(0);
                }
                spectrumIdentificationItemProteinViewTable.scrollRowToVisible(0);
                ProteinDetectionHypothesis proteinDetectionHypothesis = mzIdentMLUnmarshaller.unmarshal(ProteinDetectionHypothesis.class, (String) proteinDetectionHypothesisTable.getModel().getValueAt(row, 0));
                DBSequence dBSequence = mzIdentMLUnmarshaller.unmarshal(DBSequence.class, proteinDetectionHypothesis.getDBSequenceRef());
                List<PeptideHypothesis> peptideHypothesisList = proteinDetectionHypothesis.getPeptideHypothesis();
                String proteinSequence = "";
                if (dBSequence != null) {
                    List<CvParam> cvParamListDBSequence = dBSequence.getCvParam();
                    String scientific_name = null;
                    String protein_description = null;
                    for (int j = 0; j < cvParamListDBSequence.size(); j++) {
                        CvParam cvParam = cvParamListDBSequence.get(j);
                        if (cvParam.getName().equals("taxonomy: scientific name")) {
                            scientific_name = cvParam.getValue();
                        }
                        if (cvParam.getName().equals("protein description")) {
                            protein_description = cvParam.getValue();
                        }
                    }
                    jScientificNameValueLabel.setText(scientific_name);
                    jProteinDescriptionEditorPane.setText(protein_description);
                    proteinSequence = dBSequence.getSeq();

                    for (int i = 0; i < peptideHypothesisList.size(); i++) {
                        PeptideHypothesis peptideHypothesis = peptideHypothesisList.get(i);
                        for (int j = 0; j < peptideHypothesis.getSpectrumIdentificationItemRef().size(); j++) {
                            spectrumIdentificationItem2 = mzIdentMLUnmarshaller.unmarshal(SpectrumIdentificationItem.class, peptideHypothesis.getSpectrumIdentificationItemRef().get(j).getSpectrumIdentificationItemRef());
                        }
                        List<PeptideEvidenceRef> peptideEvidenceRefList = spectrumIdentificationItem2.getPeptideEvidenceRef();

                        for (int k = 0; k < peptideEvidenceRefList.size(); k++) {
                            PeptideEvidenceRef peptideEvidenceRef = peptideEvidenceRefList.get(k);
                            PeptideEvidence peptideEvidence = mzIdentMLUnmarshaller.unmarshal(PeptideEvidence.class, peptideEvidenceRef.getPeptideEvidenceRef());
                            if (peptideEvidence.getDBSequenceRef().equals(proteinDetectionHypothesis.getDBSequenceRef())) {
                                Peptide peptide = mzIdentMLUnmarshaller.unmarshal(Peptide.class, spectrumIdentificationItem2.getPeptideRef());
                                if (peptide != null) {
                                    List<Modification> modificationList = peptide.getModification();
                                    Modification modification;
                                    String residues = null;
                                    Integer location = null;
                                    String modificationName = null;
                                    CvParam modificationCvParam;
                                    String combine = null;
                                    if (modificationList.size() > 0) {
                                        modification = modificationList.get(0);
                                        location = modification.getLocation();
                                        if (modification.getResidues().size() > 0) {
                                            residues = modification.getResidues().get(0);
                                        }
                                        List<CvParam> modificationCvParamList = modification.getCvParam();
                                        if (modificationCvParamList.size() > 0) {
                                            modificationCvParam = modificationCvParamList.get(0);
                                            modificationName = modificationCvParam.getName();
                                        }
                                    }
                                    if (modificationName != null) {
                                        combine = modificationName;
                                    }
                                    if (residues != null) {
                                        combine = combine + " on residues: " + residues;
                                    }
                                    if (location != null) {
                                        combine = combine + " at location: " + location;
                                    }
                                    ((DefaultTableModel) spectrumIdentificationItemProteinViewTable.getModel()).addRow(new String[]{
                                                peptide.getPeptideSequence(), spectrumIdentificationItem2.getId(), combine});
                                    String find = peptide.getPeptideSequence();
                                    String replace = "<FONT COLOR=\"red\">" + find + "</FONT>";
                                    Pattern pattern = Pattern.compile(find);
                                    if (proteinSequence != null) {
                                        Matcher matcher = pattern.matcher(proteinSequence);
                                        proteinSequence = matcher.replaceAll(replace);
                                    }
                                    List<CvParam> cvParamListSpectrumIdentificationItem = spectrumIdentificationItem2.getCvParam();

                                    for (int s = 0; s < cvParamListSpectrumIdentificationItem.size(); s++) {
                                        CvParam cvParam = cvParamListSpectrumIdentificationItem.get(s);
                                        String newCol = cvParam.getName();
                                        if (newCol.equals("mascot:score")) {
                                            ((DefaultTableModel) spectrumIdentificationItemProteinViewTable.getModel()).setValueAt(Double.valueOf(cvParam.getValue()), ((DefaultTableModel) spectrumIdentificationItemProteinViewTable.getModel()).getRowCount() - 1, 3);
                                        }
                                        String accession = cvParam.getAccession();
                                        if (accession.equals("MS:1001330")
                                                || accession.equals("MS:1001172")
                                                || accession.equals("MS:1001159")
                                                || accession.equals("MS:1001328")) {
                                            ((DefaultTableModel) spectrumIdentificationItemProteinViewTable.getModel()).setValueAt(Double.valueOf(cvParam.getValue()), ((DefaultTableModel) spectrumIdentificationItemProteinViewTable.getModel()).getRowCount() - 1, 4);
                                        }


                                    }
                                }

                            }

                        }



                    }


                }
                if (proteinSequence != null) {
                    StringBuilder sb = new StringBuilder(proteinSequence);
                    StringBuilder sb_new = new StringBuilder();
                    int i = 0;

                    for (int j = 0; j < sb.length(); j++) {

                        if (i % 60 == 0 && i != 0) {
                            sb_new.append("<BR>");
                        }
                        i = i + 1;
                        sb_new.append(sb.charAt(j));
                        if (sb.charAt(j) == '<') {
                            if (sb.charAt(j + 1) == '/') {
                                sb_new.append(sb.charAt(j + 1));
                                for (int z = j + 2; z <= j + 6; z++) {
                                    sb_new.append(sb.charAt(z));
                                }
                                j = j + 6;
                            } else {
                                for (int z = j + 1; z <= j + 17; z++) {
                                    sb_new.append(sb.charAt(z));
                                }
                                j = j + 17;
                            }
                        }

                    }
                    jProteinSequenceTextPane.setText(sb_new.toString());

                }


                SimpleAttributeSet sa = new SimpleAttributeSet();
                StyleConstants.setAlignment(sa, StyleConstants.ALIGN_JUSTIFIED);

                jProteinSequenceTextPane.getStyledDocument().setParagraphAttributes(0, 60, sa, false);
            } catch (JAXBException ex) {
                Logger.getLogger(MzIdentMLViewer.class.getName()).log(Level.SEVERE, null, ex);
            }





        }
        proteinDetectionHypothesisTable.setSortable(true);
        this.setCursor(new java.awt.Cursor(java.awt.Cursor.DEFAULT_CURSOR));
    }

    /**
     * protein Ambiguity Group Table Mouse Clicked
     */
    private void proteinAmbiguityGroupTableMouseClicked(MouseEvent evt) {
        proteinDetectionHypothesisTable.setSortable(false);
        jProteinSequenceTextPane.setText("");
        jProteinSequenceTextPane.setText("");
        jScientificNameValueLabel.setText("");
        this.setCursor(new java.awt.Cursor(java.awt.Cursor.WAIT_CURSOR));
        int row = proteinAmbiguityGroupTable.getSelectedRow();



        if (row != -1) {
            try {
                while (proteinDetectionHypothesisTable.getRowCount() > 0) {
                    ((DefaultTableModel) proteinDetectionHypothesisTable.getModel()).removeRow(0);
                }
                while (((DefaultTableModel) spectrumIdentificationItemProteinViewTable.getModel()).getRowCount() > 0) {
                    ((DefaultTableModel) spectrumIdentificationItemProteinViewTable.getModel()).removeRow(0);
                }
                proteinDetectionHypothesisTable.scrollRowToVisible(0);
                String pag_id = (String) proteinAmbiguityGroupTable.getModel().getValueAt(row, 0);


                ProteinAmbiguityGroup proteinAmbiguityGroup = mzIdentMLUnmarshaller.unmarshal(ProteinAmbiguityGroup.class, pag_id);

                List<ProteinDetectionHypothesis> proteinDetectionHypothesisList = proteinAmbiguityGroup.getProteinDetectionHypothesis();
                if (proteinDetectionHypothesisList.size() > 0) {
                    for (int i = 0; i < proteinDetectionHypothesisList.size(); i++) {

                        ProteinDetectionHypothesis proteinDetectionHypothesis = proteinDetectionHypothesisList.get(i);
                        DBSequence dBSequence = mzIdentMLUnmarshaller.unmarshal(DBSequence.class, proteinDetectionHypothesis.getDBSequenceRef());
                        boolean isDecoy = checkIfProteinDetectionHypothesisIsDecoy(proteinDetectionHypothesis);
                        List<CvParam> cvParamList = proteinDetectionHypothesis.getCvParam();
                        String score = " ";
                        String number_peptide = " ";
                        for (int j = 0; j < cvParamList.size(); j++) {
                            CvParam cvParam = cvParamList.get(j);
                            if (cvParam.getName().contains("score")) {
                                score = cvParam.getValue();
                            }


                        }
                        String dBSequenceAccession = "";
                        if (dBSequence != null) {
                            dBSequenceAccession = dBSequence.getAccession();
                        }
                        if (proteinDetectionHypothesis.getPeptideHypothesis() != null) {
                            number_peptide = String.valueOf(proteinDetectionHypothesis.getPeptideHypothesis().size());
                        }
                        ((DefaultTableModel) proteinDetectionHypothesisTable.getModel()).addRow(new Object[]{
                                    proteinDetectionHypothesis.getId(),
                                    dBSequenceAccession,
                                    score,
                                    "",
                                    Integer.valueOf(number_peptide),
                                    isDecoy
                                });
                    }
                }
            } catch (JAXBException ex) {
                Logger.getLogger(MzIdentMLViewer.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
        proteinDetectionHypothesisTable.setSortable(true);
        this.setCursor(new java.awt.Cursor(java.awt.Cursor.DEFAULT_CURSOR));
    }

    private boolean checkIfProteinDetectionHypothesisIsDecoy(ProteinDetectionHypothesis proteinDetectionHypothesis) {
        boolean result = false;
        List<PeptideHypothesis> PeptideHyposthesisList = proteinDetectionHypothesis.getPeptideHypothesis();
        for (int i = 0; i < PeptideHyposthesisList.size(); i++) {
            try {
                PeptideHypothesis peptideHypothesis = PeptideHyposthesisList.get(i);
                String peptidRef = peptideHypothesis.getPeptideEvidenceRef();
                PeptideEvidence peptiedEvidence = mzIdentMLUnmarshaller.unmarshal(PeptideEvidence.class, peptidRef);
                if (peptiedEvidence.isIsDecoy()) {
                    result = true;
                    break;
                }
            } catch (JAXBException ex) {
                Logger.getLogger(MzIdentMLViewer.class.getName()).log(Level.SEVERE, null, ex);
            }

        }
        return result;
    }
    private void openMenuItemActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_openMenuItemActionPerformed
        progressBarDialog = new ProgressBarDialog(this, true);
        final Thread thread = new Thread(new Runnable() {

            @Override
            public void run() {
                progressBarDialog.setTitle("Parsing the mzid file. Please Wait...");
                progressBarDialog.setVisible(true);
            }
        }, "ProgressBarDialog");

        int returnVal = fileChooser.showOpenDialog(this);
        if (returnVal == JFileChooser.APPROVE_OPTION) {
            final File mzid_file = fileChooser.getSelectedFile();
            setTitle("MzIdentML Viewer   -  " + mzid_file.getPath());

            thread.start();
            new Thread("LoadingThread") {

                @Override
                public void run() {
                    try {
                        if (mzid_file.getPath().endsWith(".gz")) {
                            File outFile;
                            FileOutputStream fos;
                            try (GZIPInputStream gin = new GZIPInputStream(new FileInputStream(mzid_file))) {
                                outFile = new File(mzid_file.getParent(), mzid_file.getName().replaceAll("\\.gz$", ""));
                                fos = new FileOutputStream(outFile);
                                byte[] buf = new byte[100000];
                                int len;
                                while ((len = gin.read(buf)) > 0) {
                                    fos.write(buf, 0, len);
                                }
                            }
                            fos.close();

                            mzIdentMLUnmarshaller = new MzIdentMLUnmarshaller(outFile);
                        } else {
                            mzIdentMLUnmarshaller = new MzIdentMLUnmarshaller(mzid_file);
                        }




                        if (!mzIdentMLUnmarshaller.getMzIdentMLVersion().startsWith("1.1.")) {
                            progressBarDialog.setVisible(false);
                            progressBarDialog.dispose();
                            JOptionPane.showMessageDialog(null, "The file is not compatible with the Viewer: different mzIdentMl version", "mzIdentMl version", JOptionPane.INFORMATION_MESSAGE);
                            return;
                        }
                        createTables();
                        clearSummaryStats();
                        mainTabbedPane.setSelectedIndex(0);
                        firstTab = true;
                        secondTab = false;
                        thirdTab = false;
                        fourthTab = false;
                        fifthTab = false;
                        loadProteinAmbiguityGroupTable();


                        progressBarDialog.setVisible(false);
                        progressBarDialog.dispose();
                        if (proteinAmbiguityGroupTable.getRowCount() == 0) {
                            JOptionPane.showMessageDialog(null, "There is no protein view for this file", "Protein View", JOptionPane.INFORMATION_MESSAGE);
                        }
                        //loadSummaryStats();


                    } catch (OutOfMemoryError error) {
                        progressBarDialog.setVisible(false);
                        progressBarDialog.dispose();
                        Runtime.getRuntime().gc();
                        JOptionPane.showMessageDialog(null, "Out of Memory Error.", "Error", JOptionPane.ERROR_MESSAGE);

                        System.exit(0);
                    } catch (Exception ex) {
                        progressBarDialog.setVisible(false);
                        progressBarDialog.dispose();
                        System.out.println(ex.getMessage());
                    }


                }
            }.start();
        }
    }//GEN-LAST:event_openMenuItemActionPerformed
    private void loadDBSequenceTable() {
        this.setCursor(new java.awt.Cursor(java.awt.Cursor.WAIT_CURSOR));
        Iterator<DBSequence> iterDBSequence = mzIdentMLUnmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.DBSequence);
        while (iterDBSequence.hasNext()) {
            DBSequence dBSequence = iterDBSequence.next();

            String cv = "";
            if (dBSequence.getCvParam() != null) {
                if (dBSequence.getCvParam().size() > 0) {
                    cv = dBSequence.getCvParam().get(0).getValue();
                }
            }

            ((DefaultTableModel) dBSequenceTable.getModel()).addRow(new String[]{
                        dBSequence.getId(),
                        dBSequence.getAccession(),
                        dBSequence.getSearchDatabaseRef(),
                        dBSequence.getSeq(),
                        cv
                    });

        }
        dBSequenceTable.setSortable(true);
        this.setCursor(new java.awt.Cursor(java.awt.Cursor.DEFAULT_CURSOR));
    }

    private void loadProtocolData() {
        this.setCursor(new java.awt.Cursor(java.awt.Cursor.WAIT_CURSOR));
        HashMap<String, AnalysisSoftware> analysisSoftwareHashMap;
        AnalysisProtocolCollection analysisProtocolCollection;
        List<SpectrumIdentificationProtocol> spectrumIdentificationProtocol;
        ProteinDetectionProtocol proteinDetectionProtocol;

        analysisSoftwareHashMap = new HashMap();

        Iterator<AnalysisSoftware> iterAnalysisSoftware = mzIdentMLUnmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.AnalysisSoftware);
        while (iterAnalysisSoftware.hasNext()) {
            AnalysisSoftware analysisSoftware = iterAnalysisSoftware.next();
            analysisSoftwareHashMap.put(analysisSoftware.getId(), analysisSoftware);

        }

        analysisProtocolCollection = mzIdentMLUnmarshaller.unmarshal(MzIdentMLElement.AnalysisProtocolCollection);

        spectrumIdentificationProtocol = analysisProtocolCollection.getSpectrumIdentificationProtocol();
        proteinDetectionProtocol = analysisProtocolCollection.getProteinDetectionProtocol();


        if (spectrumIdentificationProtocol != null) {

            ((DefaultTableModel) protocolTable.getModel()).addRow(new String[]{
                        "<html><font color=red>Spectrum Identification Protocol</font>",
                        " "
                    });
            for (int i = 0; i < spectrumIdentificationProtocol.size(); i++) {
                SpectrumIdentificationProtocol spectrumIdentificationProtocol1 = spectrumIdentificationProtocol.get(i);
                Param param = spectrumIdentificationProtocol.get(i).getSearchType();
                if (param != null) {
                    CvParam cVParam = param.getCvParam();
                    if (cVParam != null) {
                        ((DefaultTableModel) protocolTable.getModel()).addRow(new String[]{
                                    "Type of search",
                                    cVParam.getName()
                                });
                    }
                }

                Enzymes enzymes = spectrumIdentificationProtocol.get(i).getEnzymes();
                if (enzymes != null) {
                    List<Enzyme> enzymesList = enzymes.getEnzyme();
                    for (int j = 0; j < enzymesList.size(); j++) {
                        Enzyme enzyme = enzymesList.get(j);
                        List<CvParam> cvParamList = enzyme.getEnzymeName().getCvParam();
                        for (int k = 0; k < cvParamList.size(); k++) {
                            CvParam cvParam = cvParamList.get(k);
                            ((DefaultTableModel) protocolTable.getModel()).addRow(new String[]{
                                        "Enzyme",
                                        cvParam.getName()
                                    });
                        }

                    }

                }

                ModificationParams modificationParams = spectrumIdentificationProtocol.get(i).getModificationParams();
                if (modificationParams != null) {
                    List<SearchModification> searchModificationList = modificationParams.getSearchModification();
                    for (int j = 0; j < searchModificationList.size(); j++) {


                        SearchModification searchModification = searchModificationList.get(j);
                        ((DefaultTableModel) protocolTable.getModel()).addRow(new String[]{
                                    "Search Modification",
                                    "Residues: " + searchModification.getResidues() + " Mass Delta: " + searchModification.getMassDelta()
                                });

                        if (searchModification.getCvParam() != null) {
                            if (searchModification.getCvParam().get(0) != null) {
                            }
                        }

                    }
                }


                Tolerance tolerance = spectrumIdentificationProtocol.get(i).getFragmentTolerance();
                if (tolerance != null) {
                    ((DefaultTableModel) protocolTable.getModel()).addRow(new String[]{
                                "Fragment Tolerance",
                                " "
                            });
                    List<CvParam> cvParamList = tolerance.getCvParam();
                    for (int j = 0; j < cvParamList.size(); j++) {
                        CvParam cvParam = cvParamList.get(j);
                        ((DefaultTableModel) protocolTable.getModel()).addRow(new String[]{
                                    cvParam.getName(),
                                    cvParam.getValue()
                                });
                    }
                }

                Tolerance toleranceParent = spectrumIdentificationProtocol.get(i).getParentTolerance();
                if (toleranceParent != null) {
                    ((DefaultTableModel) protocolTable.getModel()).addRow(new String[]{
                                "Parent Tolerance",
                                " "
                            });
                    List<CvParam> cvParamList = toleranceParent.getCvParam();
                    for (int j = 0; j < cvParamList.size(); j++) {
                        CvParam cvParam = cvParamList.get(j);
                        ((DefaultTableModel) protocolTable.getModel()).addRow(new String[]{
                                    cvParam.getName(),
                                    cvParam.getValue()
                                });
                    }
                }

                String sw_ref_sip = spectrumIdentificationProtocol1.getAnalysisSoftwareRef();
                if (sw_ref_sip != null && !sw_ref_sip.equals("")) {
                    ((DefaultTableModel) protocolTable.getModel()).addRow(new String[]{
                                "Analysis Software Ref",
                                sw_ref_sip
                            });
                    ((DefaultTableModel) protocolTable.getModel()).addRow(new String[]{
                                "Analysis Software name",
                                analysisSoftwareHashMap.get(sw_ref_sip).getName()
                            });
                    ((DefaultTableModel) protocolTable.getModel()).addRow(new String[]{
                                "Analysis Software version",
                                analysisSoftwareHashMap.get(sw_ref_sip).getVersion()
                            });
                    ((DefaultTableModel) protocolTable.getModel()).addRow(new String[]{
                                "Analysis Software uri",
                                analysisSoftwareHashMap.get(sw_ref_sip).getUri()
                            });

                }
            }



            ((DefaultTableModel) protocolTable.getModel()).addRow(new String[]{
                        " ",
                        " "
                    });
        }
        if (proteinDetectionProtocol != null) {

            ((DefaultTableModel) protocolTable.getModel()).addRow(new String[]{
                        "<html><font color=red>Protein Detection Protocol</font>",
                        " "
                    });
            String sw_ref_pdp = proteinDetectionProtocol.getAnalysisSoftwareRef();
            if (sw_ref_pdp != null && !sw_ref_pdp.equals("")) {
                ((DefaultTableModel) protocolTable.getModel()).addRow(new String[]{
                            "Analysis Software Ref",
                            sw_ref_pdp
                        });
                ((DefaultTableModel) protocolTable.getModel()).addRow(new String[]{
                            "Analysis Software name",
                            analysisSoftwareHashMap.get(sw_ref_pdp).getName()
                        });
                ((DefaultTableModel) protocolTable.getModel()).addRow(new String[]{
                            "Analysis Software version",
                            analysisSoftwareHashMap.get(sw_ref_pdp).getVersion()
                        });
                ((DefaultTableModel) protocolTable.getModel()).addRow(new String[]{
                            "Analysis Software uri",
                            analysisSoftwareHashMap.get(sw_ref_pdp).getUri()
                        });
            }

        }
        this.setCursor(new java.awt.Cursor(java.awt.Cursor.DEFAULT_CURSOR));
    }

    private void clearSummaryStats() {

        totalSIRLabelValue.setText("0");
        totalSIIaboveThresholdLabelValue.setText("0");
        totalSIIbelowThresholdLabelValue.setText("0");
        totalSIIaboveThresholdRankOneLabelValue.setText("0");
        percentIdentifiedSpectraLabelValue.setText("0");
        totalPeptidesaboveThresholdLabelValue.setText("0");
        totalPAGsLabelValue.setText("0");
        totalPDHsaboveThresholdLabelValue.setText("0");
        isDecoySiiValue.setText("0");
        isDecoySiiFalseValue.setText("0");


        fpSiiValue.setText("0");
        tpSiiValue.setText("0");
        fdrSiiValue.setText("0");
        siiComboBox.removeAllItems();

        while (fdrPanel.getComponents().length > 0) {
            fdrPanel.remove(0);
        }
        fdrPanel.validate();
        fdrPanel.repaint();

        while (tpEvaluePanel.getComponents().length > 0) {
            tpEvaluePanel.remove(0);
        }
        tpEvaluePanel.validate();
        tpEvaluePanel.repaint();

        while (tpQvaluePanel.getComponents().length > 0) {
            tpQvaluePanel.remove(0);
        }
        tpQvaluePanel.validate();
        tpQvaluePanel.repaint();

    }

    private void loadSummaryStats() {

        Iterator<SpectrumIdentificationList> iterspectrumIdentificationList = mzIdentMLUnmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.SpectrumIdentificationList);
        while (iterspectrumIdentificationList.hasNext()) {
            SpectrumIdentificationList spectrumIdentificationList1 = iterspectrumIdentificationList.next();

            siiComboBox.addItem(spectrumIdentificationList1.getId());
        }
        loadSpectrumIdentificationList();

    }

    private void loadSpectrumIdentificationList() {
        this.setCursor(new java.awt.Cursor(java.awt.Cursor.WAIT_CURSOR));
        try {
            sIIListPassThreshold.clear();
            sIIListBelowThreshold.clear();
            sIIListPassThresholdRankOne.clear();
            peptideListNonReduntant.clear();
            sIIListIsDecoyFalse.clear();
            sIIListIsDecoyTrue.clear();

            List<SpectrumIdentificationResult> sirListTemp = mzIdentMLUnmarshaller.unmarshal(SpectrumIdentificationList.class, siiComboBox.getSelectedItem().toString()).getSpectrumIdentificationResult();
            List<SpectrumIdentificationItem> siiListTemp = new ArrayList();
            for (int i = 0; i < sirListTemp.size(); i++) {
                SpectrumIdentificationResult spectrumIdentificationResult = sirListTemp.get(i);

                for (int j = 0; j < spectrumIdentificationResult.getSpectrumIdentificationItem().size(); j++) {
                    siiListTemp.add(spectrumIdentificationResult.getSpectrumIdentificationItem().get(j));

                }
            }


            for (int j = 0; j < siiListTemp.size(); j++) {
                SpectrumIdentificationItem spectrumIdentificationItem = siiListTemp.get(j);

                if (spectrumIdentificationItem.isPassThreshold()) {
                    sIIListPassThreshold.add(spectrumIdentificationItem);
                    String p_ref = spectrumIdentificationItem.getPeptideRef();

                    Peptide peptide = mzIdentMLUnmarshaller.unmarshal(Peptide.class, p_ref);
                    if (!peptideListNonReduntant.contains(peptide)) {
                        peptideListNonReduntant.add(peptide);
                    }
                } else {
                    sIIListBelowThreshold.add(spectrumIdentificationItem);
                }

                if (spectrumIdentificationItem.isPassThreshold() && spectrumIdentificationItem.getRank() == 1) {
                    sIIListPassThresholdRankOne.add(spectrumIdentificationItem);
                }

                boolean isdecoy = true;
                List<PeptideEvidenceRef> peRefLst = spectrumIdentificationItem.getPeptideEvidenceRef();
                for (int k = 0; k < peRefLst.size(); k++) {
                    PeptideEvidenceRef peptideEvidenceRef = peRefLst.get(k);
                    PeptideEvidence peptideEvidence = mzIdentMLUnmarshaller.unmarshal(PeptideEvidence.class, peptideEvidenceRef.getPeptideEvidenceRef());
                    if (!peptideEvidence.isIsDecoy()) {
                        sIIListIsDecoyFalse.add(spectrumIdentificationItem);
                        isdecoy = false;
                        break;
                    }
                }
                if (isdecoy) {
                    sIIListIsDecoyTrue.add(spectrumIdentificationItem);
                }
            }

            int sirSize = mzIdentMLUnmarshaller.getObjectCountForXpath(MzIdentMLElement.SpectrumIdentificationResult.getXpath());
            totalSIRLabelValue.setText("" + sirSize);

            if (sIIListPassThreshold != null) {
                totalSIIaboveThresholdLabelValue.setText(String.valueOf(sIIListPassThreshold.size()));
            }
            if (sIIListBelowThreshold != null) {
                totalSIIbelowThresholdLabelValue.setText(String.valueOf(sIIListBelowThreshold.size()));
            }
            if (sIIListPassThresholdRankOne != null) {
                totalSIIaboveThresholdRankOneLabelValue.setText(String.valueOf(sIIListPassThresholdRankOne.size()));
            }
            if (sIIListPassThresholdRankOne != null && sirSize > 0) {
                double percent = roundTwoDecimals((float) sIIListPassThresholdRankOne.size() * 100 / sirSize);
                percentIdentifiedSpectraLabelValue.setText(String.valueOf(percent) + "%");
            }
            if (peptideListNonReduntant != null) {
                totalPeptidesaboveThresholdLabelValue.setText(String.valueOf(peptideListNonReduntant.size()));
            }
            int pagSize = mzIdentMLUnmarshaller.getObjectCountForXpath(MzIdentMLElement.ProteinAmbiguityGroup.getXpath());

            totalPAGsLabelValue.setText("" + pagSize);

            if (pDHListPassThreshold != null) {
                totalPDHsaboveThresholdLabelValue.setText(String.valueOf(pDHListPassThreshold.size()));
            }
            if (sIIListIsDecoyTrue != null) {
                isDecoySiiValue.setText(String.valueOf(sIIListIsDecoyTrue.size()));
            }
            if (sIIListIsDecoyFalse != null) {
                isDecoySiiFalseValue.setText(String.valueOf(sIIListIsDecoyFalse.size()));
            }
            if (sIIListIsDecoyTrue != null) {
                falsePositiveSii = roundThreeDecimals(sIIListIsDecoyTrue.size());
                fpSiiValue.setText(String.valueOf(falsePositiveSii));
            }
            if (sIIListIsDecoyTrue != null && sIIListIsDecoyFalse != null) {
                truePositiveSii = roundThreeDecimals(sIIListPassThreshold.size() - sIIListIsDecoyTrue.size());
                tpSiiValue.setText(String.valueOf(truePositiveSii));
            }
            if (sIIListIsDecoyTrue != null && sIIListIsDecoyFalse != null) {
                if (falsePositiveSii + truePositiveSii == 0) {
                    fdrSiiValue.setText("0.0");
                } else {
                    fdrSii = roundThreeDecimals((falsePositiveSii) / (falsePositiveSii + truePositiveSii));
                    fdrSiiValue.setText(String.valueOf(fdrSii));
                }
            }
            //        thresholdValue.setText(spectrumIdentificationProtocol.get(0).getThreshold().getCvParam().get(0).getValue());
        } catch (JAXBException ex) {
            Logger.getLogger(MzIdentMLViewer.class.getName()).log(Level.SEVERE, null, ex);
        }
        this.setCursor(new java.awt.Cursor(java.awt.Cursor.DEFAULT_CURSOR));
    }

    private void exportProteinAmbiguityGrouMenuItemActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_exportProteinAmbiguityGrouMenuItemActionPerformed
        export(proteinAmbiguityGroupTable, "Export Protein Ambiguity Group");
    }//GEN-LAST:event_exportProteinAmbiguityGrouMenuItemActionPerformed

    private void exportProteinDetectionHypothesisMenuItemActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_exportProteinDetectionHypothesisMenuItemActionPerformed
        export(proteinDetectionHypothesisTable, "Export Protein Detection Hypothesis");
    }//GEN-LAST:event_exportProteinDetectionHypothesisMenuItemActionPerformed

    private void exportPeptideSpectrumMatchesMenuItemActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_exportPeptideSpectrumMatchesMenuItemActionPerformed
        export(spectrumIdentificationItemProteinViewTable, "Export Peptide Spectrum Matches");
    }//GEN-LAST:event_exportPeptideSpectrumMatchesMenuItemActionPerformed

    private void exportSpectrumIdentificationResultMenuItemActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_exportSpectrumIdentificationResultMenuItemActionPerformed
        export(spectrumIdentificationResultTable, "Export Spectrum Identification Result");
    }//GEN-LAST:event_exportSpectrumIdentificationResultMenuItemActionPerformed

    private void exportFragmentationMenuItemActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_exportFragmentationMenuItemActionPerformed
        export(fragmentationTable, "Export Fragmentation ");
    }//GEN-LAST:event_exportFragmentationMenuItemActionPerformed

    private void bCheckBoxActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_bCheckBoxActionPerformed
          updateGraph();

    }//GEN-LAST:event_bCheckBoxActionPerformed

    private void yCheckBoxActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_yCheckBoxActionPerformed
        updateGraph();
    }//GEN-LAST:event_yCheckBoxActionPerformed

    private void oneCheckBoxActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_oneCheckBoxActionPerformed
           updateGraph();
    }//GEN-LAST:event_oneCheckBoxActionPerformed

    private void twoCheckBoxActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_twoCheckBoxActionPerformed
        
           updateGraph();
    }//GEN-LAST:event_twoCheckBoxActionPerformed

    private void fdrPlotActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_fdrPlotActionPerformed
    }//GEN-LAST:event_fdrPlotActionPerformed

    private void lookAndFeelDefaultMenuItemActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_lookAndFeelDefaultMenuItemActionPerformed
        try {
            UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
        } catch (ClassNotFoundException | InstantiationException | IllegalAccessException | UnsupportedLookAndFeelException ex) {
//           
        }

        lookAndFeelDefaultMenuItem.setSelected(true);
        lookAndFeelMetalMenuItem.setSelected(false);
        lookAndFeeMotifMenuItem.setSelected(false);
        repaint();

    }//GEN-LAST:event_lookAndFeelDefaultMenuItemActionPerformed

    private void lookAndFeelMetalMenuItemActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_lookAndFeelMetalMenuItemActionPerformed
        try {
            UIManager.setLookAndFeel("javax.swing.plaf.metal.MetalLookAndFeel");
        } catch (ClassNotFoundException | InstantiationException | IllegalAccessException | UnsupportedLookAndFeelException ex) {
        }
        lookAndFeelDefaultMenuItem.setSelected(false);
        lookAndFeelMetalMenuItem.setSelected(true);
        lookAndFeeMotifMenuItem.setSelected(false);
        repaint();
    }//GEN-LAST:event_lookAndFeelMetalMenuItemActionPerformed

    private void lookAndFeeMotifMenuItemActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_lookAndFeeMotifMenuItemActionPerformed

        try {
            UIManager.setLookAndFeel("com.sun.java.swing.plaf.motif.MotifLookAndFeel");
        } catch (ClassNotFoundException | InstantiationException | IllegalAccessException | UnsupportedLookAndFeelException ex) {
        }
        lookAndFeelDefaultMenuItem.setSelected(false);
        lookAndFeelMetalMenuItem.setSelected(false);
        lookAndFeeMotifMenuItem.setSelected(true);
        repaint();
    }//GEN-LAST:event_lookAndFeeMotifMenuItemActionPerformed

    private void font12ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_font12ActionPerformed
        // TODO add your handling code here:
        changeFontSize(12);
        font12.setSelected(true);
        font14.setSelected(false);
        font16.setSelected(false);
    }//GEN-LAST:event_font12ActionPerformed

    private void font14ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_font14ActionPerformed
        // TODO add your handling code here:
        changeFontSize(14);
        font12.setSelected(false);
        font14.setSelected(true);
        font16.setSelected(false);
    }//GEN-LAST:event_font14ActionPerformed

    private void font16ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_font16ActionPerformed
        // TODO add your handling code here:
        changeFontSize(16);
        font12.setSelected(false);
        font14.setSelected(false);
        font16.setSelected(true);
    }//GEN-LAST:event_font16ActionPerformed

    private void manualCalculateActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_manualCalculateActionPerformed
        try {
            if (!getTitle().endsWith("*")) {
                setTitle(getTitle() + " *");
            }
            sIIListIsDecoyTrue.clear();
            sIIListIsDecoyFalse.clear();

            List<SpectrumIdentificationResult> sirListTemp = mzIdentMLUnmarshaller.unmarshal(SpectrumIdentificationList.class, siiComboBox.getSelectedItem().toString()).getSpectrumIdentificationResult();
            List<SpectrumIdentificationItem> siiListTemp = new ArrayList();
            for (int i = 0; i < sirListTemp.size(); i++) {
                SpectrumIdentificationResult spectrumIdentificationResult = sirListTemp.get(i);

                for (int j = 0; j < spectrumIdentificationResult.getSpectrumIdentificationItem().size(); j++) {
                    siiListTemp.add(spectrumIdentificationResult.getSpectrumIdentificationItem().get(j));

                }
            }

            boolean isdecoy = true;
            for (int j = 0; j < siiListTemp.size(); j++) {
                SpectrumIdentificationItem spectrumIdentificationItem = siiListTemp.get(j);
                List<PeptideEvidenceRef> peptideEvidenceRefList = spectrumIdentificationItem.getPeptideEvidenceRef();

                for (int k = 0; k < peptideEvidenceRefList.size(); k++) {
                    PeptideEvidenceRef peptideEvidenceRef = peptideEvidenceRefList.get(k);

                    PeptideEvidence peptideEvidence1 = mzIdentMLUnmarshaller.unmarshal(PeptideEvidence.class, peptideEvidenceRef.getPeptideEvidenceRef());
                    DBSequence dbSeq = mzIdentMLUnmarshaller.unmarshal(DBSequence.class, peptideEvidence1.getDBSequenceRef());

                    if (!dbSeq.getAccession().startsWith(manualDecoyPrefixValue.getText())) {
                        sIIListIsDecoyFalse.add(spectrumIdentificationItem);
                        isdecoy = false;
                        break;

                    }

                }
                if (isdecoy) {
                    sIIListIsDecoyTrue.add(spectrumIdentificationItem);
                }


            }




            if (sIIListIsDecoyTrue != null) {
                isDecoySiiValue.setText(String.valueOf(sIIListIsDecoyTrue.size()));
            }
            if (sIIListIsDecoyFalse != null) {
                isDecoySiiFalseValue.setText(String.valueOf(sIIListIsDecoyFalse.size()));
            }
            if (sIIListIsDecoyTrue != null) {
                falsePositiveSii = roundThreeDecimals(sIIListIsDecoyTrue.size() / Double.valueOf(manualDecoyRatioValue.getText().trim()));
                fpSiiValue.setText(String.valueOf(falsePositiveSii));
            }
            if (sIIListIsDecoyTrue != null && sIIListIsDecoyFalse != null) {
                truePositiveSii = roundThreeDecimals(sIIListPassThreshold.size() - sIIListIsDecoyTrue.size());
                tpSiiValue.setText(String.valueOf(truePositiveSii));
            }
            if (sIIListIsDecoyTrue != null && sIIListIsDecoyFalse != null) {
                fdrSii = roundThreeDecimals(falsePositiveSii / (falsePositiveSii + truePositiveSii));
                fdrSiiValue.setText(String.valueOf(fdrSii));
            }

            while (fdrPanel.getComponents().length > 0) {
                fdrPanel.remove(0);
            }
            falseDiscoveryRate = new FalseDiscoveryRate(mzIdentMLUnmarshaller, manualDecoyRatioValue.getText(), manualDecoyPrefixValue.getText(), siiComboBox.getSelectedIndex());
            falseDiscoveryRate.computeFDRusingJonesMethod();



            // FDR Graph
            XYSeriesCollection datasetFDR = new XYSeriesCollection();
            final XYSeries dataFDR = new XYSeries("FDR", false);

            for (int i = 0; i < falseDiscoveryRate.getSorted_evalues().size(); i++) {
                dataFDR.add(Math.log10(falseDiscoveryRate.getSorted_evalues().get(i)), falseDiscoveryRate.getSorted_estimatedFDR().get(i));

            }

            final XYSeries dataFDRQvalue = new XYSeries("Q-value", false);

            for (int i = 0; i < falseDiscoveryRate.getSorted_evalues().size(); i++) {
                dataFDRQvalue.add(Math.log10(falseDiscoveryRate.getSorted_evalues().get(i)), falseDiscoveryRate.getSorted_qValues().get(i));

            }

            final XYSeries dataFDRSimple = new XYSeries("Simple FDR", false);

            for (int i = 0; i < falseDiscoveryRate.getSorted_evalues().size(); i++) {
                dataFDRSimple.add(Math.log10(falseDiscoveryRate.getSorted_evalues().get(i)), falseDiscoveryRate.getSorted_simpleFDR().get(i));

            }




            datasetFDR.addSeries(dataFDR);
            datasetFDR.addSeries(dataFDRQvalue);
            datasetFDR.addSeries(dataFDRSimple);
            final JFreeChart chartFDR = createFDRChart(datasetFDR);
            final ChartPanel chartPanelFDR = new ChartPanel(chartFDR);

            chartPanelFDR.setPreferredSize(new java.awt.Dimension(257, 255));
            fdrPanel.add(chartPanelFDR);
            JScrollPane jFDRPane = new JScrollPane(chartPanelFDR);
            fdrPanel.setLayout(new java.awt.BorderLayout());
            fdrPanel.add(jFDRPane);

            // TP vs Evalue Graph
            XYSeriesCollection datasetTpQvalue = new XYSeriesCollection();
            final XYSeries dataTpQvalue = new XYSeries("TP", false);

            for (int i = 0; i < falseDiscoveryRate.getSorted_evalues().size(); i++) {
                dataTpQvalue.add(Math.log10(falseDiscoveryRate.getSorted_evalues().get(i)), falseDiscoveryRate.getTP().get(i));

            }

            final XYSeries dataFpQvalue = new XYSeries("FP", false);

            for (int i = 0; i < falseDiscoveryRate.getSorted_evalues().size(); i++) {
                dataFpQvalue.add(Math.log10(falseDiscoveryRate.getSorted_evalues().get(i)), falseDiscoveryRate.getFP().get(i));

            }


            datasetTpQvalue.addSeries(dataTpQvalue);
            datasetTpQvalue.addSeries(dataFpQvalue);
            final JFreeChart chartTpQvalue = createTpQvalueChart(datasetTpQvalue);
            final ChartPanel chartPanelTpQvalue = new ChartPanel(chartTpQvalue);

            chartPanelTpQvalue.setPreferredSize(new java.awt.Dimension(257, 255));
            tpEvaluePanel.add(chartPanelTpQvalue);
            JScrollPane jTpQvaluePane = new JScrollPane(chartPanelTpQvalue);
            tpEvaluePanel.setLayout(new java.awt.BorderLayout());
            tpEvaluePanel.add(jTpQvaluePane);


            // TP vs Qvalue
            XYSeriesCollection datasetTpQvalueCollection = new XYSeriesCollection();
            final XYSeries dataTpQValueSeries = new XYSeries("data", false);

            for (int i = 0; i < falseDiscoveryRate.getSorted_evalues().size(); i++) {
                dataTpQValueSeries.add(falseDiscoveryRate.getSorted_qValues().get(i), falseDiscoveryRate.getTP().get(i));

            }


            datasetTpQvalueCollection.addSeries(dataTpQValueSeries);
            final JFreeChart chartTpQvalueChart = createTpQvalue(datasetTpQvalueCollection);
            final ChartPanel chartPanelTpSimpleFDR = new ChartPanel(chartTpQvalueChart);

            chartPanelTpSimpleFDR.setPreferredSize(new java.awt.Dimension(257, 255));
            tpQvaluePanel.add(chartPanelTpSimpleFDR);
            JScrollPane jTpSimpleFDRPane = new JScrollPane(chartPanelTpSimpleFDR);
            tpQvaluePanel.setLayout(new java.awt.BorderLayout());
            tpQvaluePanel.add(jTpSimpleFDRPane);

            repaint();
        } catch (JAXBException ex) {
            Logger.getLogger(MzIdentMLViewer.class.getName()).log(Level.SEVERE, null, ex);
        }
    }//GEN-LAST:event_manualCalculateActionPerformed

    private void manualDecoyActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_manualDecoyActionPerformed
        if (manualDecoy.isSelected()) {
            manualDecoyPrefixValue.setEnabled(true);
            manualDecoyRatioValue.setEnabled(true);
            //  manualCalculate.setEnabled(true);
        } else {
            manualDecoyPrefixValue.setEnabled(false);
            manualDecoyRatioValue.setEnabled(false);
            // manualCalculate.setEnabled(false);
        }
}//GEN-LAST:event_manualDecoyActionPerformed

    private void saveMenuItemActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_saveMenuItemActionPerformed

        JFileChooser chooser = new JFileChooser();
        chooser.setFileFilter(new MzIdentMLFilter());
        chooser.setMultiSelectionEnabled(false);
        chooser.setDialogTitle("Save");

        File selectedFile;

        int returnVal = chooser.showSaveDialog(this);

        if (returnVal == JFileChooser.APPROVE_OPTION) {

            selectedFile = chooser.getSelectedFile();

            if (!selectedFile.getName().toLowerCase().endsWith(".mzid")) {
                selectedFile = new File(selectedFile.getAbsolutePath() + ".mzid");
            }

            while (selectedFile.exists()) {
                int option = JOptionPane.showConfirmDialog(this,
                        "The  file " + chooser.getSelectedFile().getName()
                        + " already exists. Replace file?",
                        "Replace File?", JOptionPane.YES_NO_CANCEL_OPTION);

                if (option == JOptionPane.NO_OPTION) {
                    chooser = new JFileChooser();
                    chooser.setFileFilter(new CsvFileFilter());
                    chooser.setMultiSelectionEnabled(false);
                    chooser.setDialogTitle("Save");

                    returnVal = chooser.showSaveDialog(this);

                    if (returnVal == JFileChooser.CANCEL_OPTION) {
                        return;
                    } else {
                        selectedFile = chooser.getSelectedFile();

                        if (!selectedFile.getName().toLowerCase().endsWith(".mzid")) {
                            selectedFile = new File(selectedFile.getAbsolutePath() + ".mzid");
                        }
                    }
                } else { // YES option
                    break;
                }
            }

            this.setCursor(new java.awt.Cursor(java.awt.Cursor.WAIT_CURSOR));


            try {

                selectedFile = chooser.getSelectedFile();

                if (!selectedFile.getName().toLowerCase().endsWith(".mzid")) {
                    selectedFile = new File(selectedFile.getAbsolutePath() + ".mzid");
                }

                if (selectedFile.exists()) {
                    selectedFile.delete();
                }

                selectedFile.createNewFile();
                try (FileWriter f = new FileWriter(selectedFile)) {
                    if (falseDiscoveryRate != null) {
                        falseDiscoveryRate.writeToMzIdentMLFile(selectedFile.getPath());
                    }
                }



            } catch (IOException ex) {
                JOptionPane.showMessageDialog(this,
                        "An error occured when exporting the spectra file details.",
                        "Error Exporting",
                        JOptionPane.ERROR_MESSAGE);
            }

            this.setCursor(new java.awt.Cursor(java.awt.Cursor.DEFAULT_CURSOR));
        }

    }//GEN-LAST:event_saveMenuItemActionPerformed

    private void siiComboBoxActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_siiComboBoxActionPerformed
        if (siiComboBox.getSelectedIndex() != -1) {
            loadSpectrumIdentificationList();
        }
    }//GEN-LAST:event_siiComboBoxActionPerformed

    private void exportFDRActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_exportFDRActionPerformed
        exportFDR();
    }//GEN-LAST:event_exportFDRActionPerformed

    private void mainTabbedPaneMouseClicked(java.awt.event.MouseEvent evt) {//GEN-FIRST:event_mainTabbedPaneMouseClicked
        // compare which tab is selected & boolean
        if (mainTabbedPane.getSelectedIndex() == 1 && !secondTab && mzIdentMLUnmarshaller != null) {
            loadSpectrumIdentificationResultTable();
            secondTab = true;
        }

        if (mainTabbedPane.getSelectedIndex() == 2 && !thirdTab && mzIdentMLUnmarshaller != null) {
            loadDBSequenceTable();
            thirdTab = true;
        }

        if (mainTabbedPane.getSelectedIndex() == 3 && !fourthTab && mzIdentMLUnmarshaller != null) {
            loadSummaryStats();
            fourthTab = true;
        }

        if (mainTabbedPane.getSelectedIndex() == 4 && !fifthTab && mzIdentMLUnmarshaller != null) {
            loadProtocolData();
            fifthTab = true;
        }

    }//GEN-LAST:event_mainTabbedPaneMouseClicked
    private void changeFontSize(int i) {
        proteinAmbiguityGroupTable.setFont(new Font("Serif", Font.PLAIN, i));
        proteinDetectionHypothesisTable.setFont(new Font("Serif", Font.PLAIN, i));
        spectrumIdentificationItemProteinViewTable.setFont(new Font("Serif", Font.PLAIN, i));
        spectrumIdentificationResultTable.setFont(new Font("Serif", Font.PLAIN, i));
        spectrumIdentificationItemTable.setFont(new Font("Serif", Font.PLAIN, i));
        fragmentationTable.setFont(new Font("Serif", Font.PLAIN, i));
        protocolTable.setFont(new Font("Serif", Font.PLAIN, i));
    }

    private JFreeChart createFDRChart(final XYDataset dataset) {
        final JFreeChart chart = ChartFactory.createScatterPlot(
                "FDR", // chart title
                "log10(e-value)", // x axis label
                "", // y axis label
                dataset, // data
                PlotOrientation.VERTICAL,
                true, // include legend
                true, // tooltips
                false // urls
                );
        XYPlot plot = (XYPlot) chart.getPlot();
        XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer();
        renderer.setSeriesLinesVisible(0, true);
        plot.setRenderer(renderer);
        return chart;
    }

    private JFreeChart createTpQvalueChart(final XYDataset dataset) {
        final JFreeChart chart = ChartFactory.createScatterPlot(
                "TP vs FP vs E-value", // chart title
                "log10(e-value)", // x axis label
                "", // y axis label
                dataset, // data
                PlotOrientation.VERTICAL,
                true, // include legend
                true, // tooltips
                false // urls
                );
        XYPlot plot = (XYPlot) chart.getPlot();
        XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer();
        renderer.setSeriesLinesVisible(0, true);
        plot.setRenderer(renderer);
        return chart;
    }

    private JFreeChart createTpQvalue(final XYDataset dataset) {
        final JFreeChart chart = ChartFactory.createScatterPlot(
                "TP vs Q-value", // chart title
                "Q-value", // x axis label
                "TP value", // y axis label
                dataset, // data
                PlotOrientation.VERTICAL,
                true, // include legend
                true, // tooltips
                false // urls
                );
        XYPlot plot = (XYPlot) chart.getPlot();
        XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer();
        renderer.setSeriesLinesVisible(0, true);
        plot.setRenderer(renderer);
        return chart;
    }

   

    private void export(JXTable table, String boxTitle) {
        JFileChooser chooser = new JFileChooser();
        chooser.setFileFilter(new CsvFileFilter());
        chooser.setMultiSelectionEnabled(false);
        chooser.setDialogTitle(boxTitle);

        File selectedFile;

        int returnVal = chooser.showSaveDialog(this);

        if (returnVal == JFileChooser.APPROVE_OPTION) {

            selectedFile = chooser.getSelectedFile();

            if (!selectedFile.getName().toLowerCase().endsWith(".csv")) {
                selectedFile = new File(selectedFile.getAbsolutePath() + ".csv");
            }

            while (selectedFile.exists()) {
                int option = JOptionPane.showConfirmDialog(this,
                        "The  file " + chooser.getSelectedFile().getName()
                        + " already exists. Replace file?",
                        "Replace File?", JOptionPane.YES_NO_CANCEL_OPTION);

                if (option == JOptionPane.NO_OPTION) {
                    chooser = new JFileChooser();
                    chooser.setFileFilter(new CsvFileFilter());
                    chooser.setMultiSelectionEnabled(false);
                    chooser.setDialogTitle(boxTitle);

                    returnVal = chooser.showSaveDialog(this);

                    if (returnVal == JFileChooser.CANCEL_OPTION) {
                        return;
                    } else {
                        selectedFile = chooser.getSelectedFile();

                        if (!selectedFile.getName().toLowerCase().endsWith(".csv")) {
                            selectedFile = new File(selectedFile.getAbsolutePath() + ".csv");
                        }
                    }
                } else { // YES option
                    break;
                }
            }

            this.setCursor(new java.awt.Cursor(java.awt.Cursor.WAIT_CURSOR));


            try {

                selectedFile = chooser.getSelectedFile();

                if (!selectedFile.getName().toLowerCase().endsWith(".csv")) {
                    selectedFile = new File(selectedFile.getAbsolutePath() + ".csv");
                }

                if (selectedFile.exists()) {
                    selectedFile.delete();
                }

                selectedFile.createNewFile();
                try (FileWriter f = new FileWriter(selectedFile)) {
                    for (int j = 0; j < table.getColumnCount() - 1; j++) {
                        f.write(table.getColumnName(j) + ",");
                    }

                    f.write(table.getColumnName(table.getColumnCount() - 1) + "\r\n");

                    // add the table contents
                    for (int i = 0; i < table.getRowCount(); i++) {
                        for (int j = 0; j < table.getColumnCount() - 1; j++) {
                            f.write(table.getValueAt(i, j) + ",");
                        }

                        f.write(table.getValueAt(i, table.getColumnCount() - 1) + "\r\n");
                    }
                }



            } catch (IOException ex) {
                JOptionPane.showMessageDialog(this,
                        "An error occured when exporting the spectra file details.",
                        "Error Exporting",
                        JOptionPane.ERROR_MESSAGE);
            }

            this.setCursor(new java.awt.Cursor(java.awt.Cursor.DEFAULT_CURSOR));
        }
    }

    private void exportFDR() {
        JFileChooser chooser = new JFileChooser();
        chooser.setFileFilter(new CsvFileFilter());
        chooser.setMultiSelectionEnabled(false);
        chooser.setDialogTitle("Export FDR");

        File selectedFile;

        int returnVal = chooser.showSaveDialog(this);

        if (returnVal == JFileChooser.APPROVE_OPTION) {

            selectedFile = chooser.getSelectedFile();

            if (!selectedFile.getName().toLowerCase().endsWith(".csv")) {
                selectedFile = new File(selectedFile.getAbsolutePath() + ".csv");
            }

            while (selectedFile.exists()) {
                int option = JOptionPane.showConfirmDialog(this,
                        "The  file " + chooser.getSelectedFile().getName()
                        + " already exists. Replace file?",
                        "Replace File?", JOptionPane.YES_NO_CANCEL_OPTION);

                if (option == JOptionPane.NO_OPTION) {
                    chooser = new JFileChooser();
                    chooser.setFileFilter(new CsvFileFilter());
                    chooser.setMultiSelectionEnabled(false);
                    chooser.setDialogTitle("Export FDR");

                    returnVal = chooser.showSaveDialog(this);

                    if (returnVal == JFileChooser.CANCEL_OPTION) {
                        return;
                    } else {
                        selectedFile = chooser.getSelectedFile();

                        if (!selectedFile.getName().toLowerCase().endsWith(".csv")) {
                            selectedFile = new File(selectedFile.getAbsolutePath() + ".csv");
                        }
                    }
                } else { // YES option
                    break;
                }
            }

            this.setCursor(new java.awt.Cursor(java.awt.Cursor.WAIT_CURSOR));


            try {

                selectedFile = chooser.getSelectedFile();

                if (!selectedFile.getName().toLowerCase().endsWith(".csv")) {
                    selectedFile = new File(selectedFile.getAbsolutePath() + ".csv");
                }

                if (selectedFile.exists()) {
                    selectedFile.delete();
                }

                selectedFile.createNewFile();
                try (FileWriter f = new FileWriter(selectedFile)) {
                    String outStrHead = "sorted_spectrumResult.get(i)\tsorted_peptideNames.get(i) \t sorted_decoyOrNot.get(i) \t  sorted_evalues.get(i).toString() \t + sorted_scores.get(i).toString() \t estimated_simpleFDR.get(i) \t estimated_qvalue.get(i) \t estimated_fdrscore.get(i) \n";
                    f.write(outStrHead);
                    for (int i = 0; i < falseDiscoveryRate.getSorted_evalues().size(); i++) {


                        String outStr = falseDiscoveryRate.getSorted_spectrumResult().get(i) + "\t"
                                + falseDiscoveryRate.getSorted_peptideNames().get(i) + "\t" + falseDiscoveryRate.getSorted_decoyOrNot().get(i) + "\t"
                                //                    + sorted_evalues.get(i).toString() + "\t" + sorted_scores.get(i).toString() + "\t"
                                + falseDiscoveryRate.getSorted_simpleFDR().get(i) + "\t" + falseDiscoveryRate.getSorted_qValues().get(i) + "\t"
                                + falseDiscoveryRate.getSorted_estimatedFDR().get(i) + "\n";

                        f.write(outStr);
                    }
                }



            } catch (IOException ex) {
                JOptionPane.showMessageDialog(this,
                        "An error occured when exporting the spectra file details.",
                        "Error Exporting",
                        JOptionPane.ERROR_MESSAGE);
            }

            this.setCursor(new java.awt.Cursor(java.awt.Cursor.DEFAULT_CURSOR));
        }
    }

    /**
     * @param args the command line arguments
     */
    public static void main(String args[]) {
        java.awt.EventQueue.invokeLater(new Runnable() {

            @Override
            public void run() {

                new MzIdentMLViewer().setVisible(true);
            }
        });
    }
    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JMenuItem aboutMenuItem;
    private javax.swing.JCheckBox bCheckBox;
    private javax.swing.JPanel dBSequencePanel;
    private javax.swing.JMenuItem exitMenuItem;
    private javax.swing.JMenuItem exportFDR;
    private javax.swing.JMenuItem exportFragmentationMenuItem;
    private javax.swing.JMenu exportMenu;
    private javax.swing.JMenuItem exportPeptideSpectrumMatchesMenuItem;
    private javax.swing.JMenuItem exportProteinAmbiguityGrouMenuItem;
    private javax.swing.JMenuItem exportProteinDetectionHypothesisMenuItem;
    private javax.swing.JPopupMenu.Separator exportSeparator1;
    private javax.swing.JPopupMenu.Separator exportSeparator2;
    private javax.swing.JMenuItem exportSpectrumIdentificationResultMenuItem;
    private javax.swing.JPanel fdrPanel;
    private javax.swing.JLabel fdrProteinsLabel;
    private javax.swing.JLabel fdrProteinsValue;
    private javax.swing.JLabel fdrSiiLabel;
    private javax.swing.JLabel fdrSiiValue;
    private javax.swing.JFileChooser fileChooser;
    private javax.swing.JMenu fileMenu;
    private javax.swing.JRadioButtonMenuItem font12;
    private javax.swing.JRadioButtonMenuItem font14;
    private javax.swing.JMenuItem font16;
    private javax.swing.JMenu fontMenu;
    private javax.swing.JLabel fpProteinsLabel;
    private javax.swing.JLabel fpProteinsValue;
    private javax.swing.JLabel fpSiiLabel;
    private javax.swing.JLabel fpSiiValue;
    private javax.swing.JPanel globalStatisticsPanel;
    private javax.swing.JMenu helpMenu;
    private javax.swing.JLabel isDecoySii;
    private javax.swing.JLabel isDecoySiiFalse;
    private javax.swing.JLabel isDecoySiiFalseValue;
    private javax.swing.JLabel isDecoySiiValue;
    private javax.swing.JPanel jExperimentalFilterPanel;
    private javax.swing.JPanel jFragmentationPanel;
    private javax.swing.JPanel jGraph;
    private javax.swing.JMenuBar jMenuBar;
    private javax.swing.JPanel jPeptideEvidencePanel;
    private javax.swing.JPanel jProteinAmbiguityGroupPanel;
    private javax.swing.JEditorPane jProteinDescriptionEditorPane;
    private javax.swing.JPanel jProteinDescriptionPanel;
    private javax.swing.JScrollPane jProteinDescriptionScrollPane;
    private javax.swing.JPanel jProteinDetectionHypothesisPanel;
    private javax.swing.JPanel jProteinInfoPanel;
    private javax.swing.JPanel jProteinSequencePanel;
    private javax.swing.JScrollPane jProteinSequenceScrollPane;
    private javax.swing.JTextPane jProteinSequenceTextPane;
    private javax.swing.JLabel jScientificNameLabel;
    private javax.swing.JLabel jScientificNameValueLabel;
    private javax.swing.JSeparator jSeparator1;
    private javax.swing.JSeparator jSeparator2;
    private javax.swing.JPopupMenu.Separator jSeparator3;
    private javax.swing.JSeparator jSeparator4;
    private javax.swing.JPanel jSpectrumIdentificationItemPanel;
    private javax.swing.JPanel jSpectrumIdentificationItemProteinPanel;
    private javax.swing.JPanel jSpectrumIdentificationResultPanel;
    private javax.swing.JPanel jSpectrumPanel;
    private javax.swing.JPanel jTheoreticalFilterPanel;
    private javax.swing.JRadioButtonMenuItem lookAndFeeMotifMenuItem;
    private javax.swing.JRadioButtonMenuItem lookAndFeelDefaultMenuItem;
    private javax.swing.JMenu lookAndFeelMenu;
    private javax.swing.JRadioButtonMenuItem lookAndFeelMetalMenuItem;
    private javax.swing.JPanel mainPanel;
    private javax.swing.JTabbedPane mainTabbedPane;
    private javax.swing.JButton manualCalculate;
    private javax.swing.JCheckBox manualDecoy;
    private javax.swing.JLabel manualDecoyLabel;
    private javax.swing.JLabel manualDecoyPrefix;
    private javax.swing.JTextField manualDecoyPrefixValue;
    private javax.swing.JLabel manualDecoyRatio;
    private javax.swing.JTextField manualDecoyRatioValue;
    private javax.swing.JCheckBox oneCheckBox;
    private javax.swing.JMenuItem openMenuItem;
    private javax.swing.JMenu optionsMenu;
    private javax.swing.JLabel percentIdentifiedSpectraLabel;
    private javax.swing.JLabel percentIdentifiedSpectraLabelValue;
    private javax.swing.JPanel proteinDBViewPanel;
    private javax.swing.JPanel proteinViewPanel;
    private javax.swing.JPanel protocolPanel;
    private javax.swing.JPanel protocolSummaryPanel;
    private javax.swing.JMenuItem saveMenuItem;
    private javax.swing.JComboBox siiComboBox;
    private javax.swing.JLabel siiListLabel;
    private javax.swing.JPanel spectrumViewPanel;
    private javax.swing.JPanel summaryPanel;
    private javax.swing.JButton thresholdButton;
    private javax.swing.JLabel thresholdLabel;
    private javax.swing.JTextField thresholdValue;
    private javax.swing.JLabel totalPAGsLabel;
    private javax.swing.JLabel totalPAGsLabelValue;
    private javax.swing.JLabel totalPDHsLabel;
    private javax.swing.JLabel totalPDHsaboveThresholdLabel;
    private javax.swing.JLabel totalPDHsaboveThresholdLabelValue;
    private javax.swing.JLabel totalPeptidesaboveThresholdLabel;
    private javax.swing.JLabel totalPeptidesaboveThresholdLabelValue;
    private javax.swing.JLabel totalSIILabel;
    private javax.swing.JLabel totalSIIaboveThresholdLabel;
    private javax.swing.JLabel totalSIIaboveThresholdLabelValue;
    private javax.swing.JLabel totalSIIaboveThresholdRankOneLabel;
    private javax.swing.JLabel totalSIIaboveThresholdRankOneLabelValue;
    private javax.swing.JLabel totalSIIbelowThresholdLabel;
    private javax.swing.JLabel totalSIIbelowThresholdLabelValue;
    private javax.swing.JLabel totalSIRLabel;
    private javax.swing.JLabel totalSIRLabelValue;
    private javax.swing.JPanel tpEvaluePanel;
    private javax.swing.JLabel tpProteinsLabel;
    private javax.swing.JLabel tpProteinsValue;
    private javax.swing.JPanel tpQvaluePanel;
    private javax.swing.JLabel tpSiiLabel;
    private javax.swing.JLabel tpSiiValue;
    private javax.swing.JCheckBox twoCheckBox;
    private javax.swing.JCheckBox yCheckBox;
    // End of variables declaration//GEN-END:variables
}
