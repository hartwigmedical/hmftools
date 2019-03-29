package com.hartwig.hmftools.svanalysis.annotators;

import static com.hartwig.hmftools.common.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.INV;
import static com.hartwig.hmftools.svanalysis.analysis.CNAnalyser.CENTROMERE_CN;
import static com.hartwig.hmftools.svanalysis.analysis.CNAnalyser.P_ARM_TELOMERE_CN;
import static com.hartwig.hmftools.svanalysis.analysis.CNAnalyser.Q_ARM_TELOMERE_CN;
import static com.hartwig.hmftools.svanalysis.analysis.ClusterAnalyser.SHORT_TI_LENGTH;
import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.CHROMOSOME_ARM_P;
import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.getChromosomalArm;
import static com.hartwig.hmftools.svanalysis.types.SvCluster.RESOLVED_TYPE_LOW_QUALITY;
import static com.hartwig.hmftools.svannotation.SvGeneTranscriptCollection.PSEUDO_GENE_DATA_EXON_LENGTH;
import static com.hartwig.hmftools.svannotation.SvGeneTranscriptCollection.PSEUDO_GENE_DATA_EXON_MAX;
import static com.hartwig.hmftools.svannotation.SvGeneTranscriptCollection.PSEUDO_GENE_DATA_EXON_RANK;
import static com.hartwig.hmftools.svannotation.SvGeneTranscriptCollection.PSEUDO_GENE_DATA_TRANS_ID;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.sql.Timestamp;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.SortedSetMultimap;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverType;
import com.hartwig.hmftools.common.genepanel.HmfGenePanelSupplier;
import com.hartwig.hmftools.common.purple.copynumber.CopyNumberMethod;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.gene.ImmutableGeneCopyNumber;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;
import com.hartwig.hmftools.common.region.HmfTranscriptRegion;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneAnnotation;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;
import com.hartwig.hmftools.svanalysis.types.DriverGeneData;
import com.hartwig.hmftools.svanalysis.types.SvBreakend;
import com.hartwig.hmftools.svanalysis.types.SvCluster;
import com.hartwig.hmftools.svanalysis.types.SvLOH;
import com.hartwig.hmftools.svanalysis.types.SvLinkedPair;
import com.hartwig.hmftools.svanalysis.types.SvVarData;
import com.hartwig.hmftools.svannotation.SvGeneTranscriptCollection;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class DriverGeneAnnotator
{
    private static final Logger LOGGER = LogManager.getLogger(DriverGeneAnnotator.class);

    final DatabaseAccess mDbAccess;
    SvGeneTranscriptCollection mGeneTranscriptCollection;

    private Map<String, HmfTranscriptRegion> mAllGenesMap;
    private List<DriverCatalog> mDriverCatalog;
    private List<GeneCopyNumber> mGeneCopyNumberData;
    private BufferedWriter mFileWriter;
    private String mOutputDir;
    private double mSamplePloidy;
    private List<DriverGeneData> mDriverGeneDataList;

    private PerformanceCounter mPerfCounter;

    // references only
    private String mSampleId;
    private Map<String, List<SvBreakend>> mChrBreakendMap;
    private Map<String, List<SvLOH>> mSampleLohMap;
    private List<SvLOH> mSampleLOHData;
    private Map<String, double[]> mChrCopyNumberMap;
    private Map<String, List<GeneCopyNumber>> mSampleGeneCopyNumberMap;
    private VisualiserWriter mVisWriter;

    // temp
    private boolean mWriteMatchedGeneCopyNumber;
    private BufferedWriter mGCNFileWriter;

    public DriverGeneAnnotator(DatabaseAccess dbAccess, SvGeneTranscriptCollection geneTranscriptCollection, final String outputDir)
    {
        mDbAccess = dbAccess;
        mGeneTranscriptCollection = geneTranscriptCollection;

        mDriverCatalog = Lists.newArrayList();
        mGeneCopyNumberData = Lists.newArrayList();
        mSampleLOHData = Lists.newArrayList();
        mDriverGeneDataList = Lists.newArrayList();
        mSampleGeneCopyNumberMap = new HashMap();
        mChrCopyNumberMap = null;
        mFileWriter = null;
        mOutputDir = outputDir;
        mSamplePloidy = 0;

        mWriteMatchedGeneCopyNumber = false;
        mGCNFileWriter = null;
        mVisWriter = null;

        mPerfCounter = new PerformanceCounter("Drivers");
    }

    private static final String WRITE_GCN_DATA = "write_gcn_data";
    private static final String GCN_DATA_FILE = "gcn_data_file";

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(WRITE_GCN_DATA, false, "Write a cache of driver-matched gene copy number data");
        options.addOption(GCN_DATA_FILE, true, "Cache of driver-matched gene copy number data");
    }

    public boolean loadConfig(final CommandLine cmd)
    {
        mWriteMatchedGeneCopyNumber = cmd.hasOption(WRITE_GCN_DATA);

        initialiseGeneData(cmd.getOptionValue(GCN_DATA_FILE,""));

        return true;
    }

    public void setLohData(final Map<String, List<SvLOH>> lohData) { mSampleLohMap = lohData; }public void writeMatchedGeneCopyNumber() { mWriteMatchedGeneCopyNumber = true; }
    public void setChrCopyNumberMap(Map<String, double[]> chrCopyNumberMap) { mChrCopyNumberMap = chrCopyNumberMap; }
    public void setSamplePloidy(double ploidy)
    {
        mSamplePloidy = ploidy;
    }

    public void setVisWriter(VisualiserWriter writer) { mVisWriter = writer; }
    public final List<DriverGeneData> getDriverGeneDataList() { return mDriverGeneDataList; }

    private void initialiseGeneData(final String geneCopyNumberFile)
    {
        SortedSetMultimap<String, HmfTranscriptRegion> genesByChromosomeMap = HmfGenePanelSupplier.allGenesPerChromosomeMap37();

        mAllGenesMap = Maps.newHashMap();
        for (final HmfTranscriptRegion region : genesByChromosomeMap.values())
        {
            mAllGenesMap.put(region.gene(), region);
        }

        loadGeneCopyNumberDataFile(geneCopyNumberFile);
    }

    private final HmfTranscriptRegion findRegion(final String geneName)
    {
        return mAllGenesMap.get(geneName);
    }

    private void loadDriverCatalog(final String sampleId)
    {
        mDriverCatalog.clear();
        mDriverGeneDataList.clear();
        mDriverCatalog.addAll(mDbAccess.readDriverCatalog(sampleId));

        LOGGER.debug("sample({}) retrieved {} driver gene records", sampleId, mDriverCatalog.size());
    }

    private void loadGeneCopyNumberData(final String sampleId)
    {
        mGeneCopyNumberData.clear();

        if(!mSampleGeneCopyNumberMap.isEmpty())
        {
            final List<GeneCopyNumber> gcnDataList = mSampleGeneCopyNumberMap.get(sampleId);

            if(gcnDataList != null)
                mGeneCopyNumberData.addAll(gcnDataList);

            return;
        }

        if(!mWriteMatchedGeneCopyNumber)
        {
            mGeneCopyNumberData = mDbAccess.readGeneCopynumbers(sampleId);
            return;
        }

        LOGGER.info("sample({}) writing gene copy number data", mSampleId);

        List<GeneCopyNumber> gcnDataList = mDbAccess.readGeneCopynumbers(sampleId);

        for(final DriverCatalog driverData : mDriverCatalog)
        {
            for(final GeneCopyNumber gcnData : gcnDataList)
            {
                if(gcnData.gene().equals(driverData.gene()))
                {
                    writeGeneCopyNumberData(gcnData);
                    mGeneCopyNumberData.add(gcnData);
                    break;
                }
            }
        }

        LOGGER.info("sample({}) drivers matched({}/{}) from {} gene copy number data",
                mSampleId, mGeneCopyNumberData.size(), mDriverCatalog.size(), gcnDataList.size());
    }

    public void annotateSVs(final String sampleId, final List<SvCluster> clusters, final Map<String, List<SvBreakend>> chrBreakendMap)
    {
        mPerfCounter.start();

        mSampleId = sampleId;
        mChrBreakendMap = chrBreakendMap;

        checkPseudoGeneAnnotations(clusters);

        loadDriverCatalog(sampleId);

        if (mDriverCatalog.isEmpty())
            return;

        loadGeneCopyNumberData(sampleId);

        mSampleLOHData.clear();
        List<SvLOH> sampleLohEvents = mSampleLohMap.get(sampleId);

        if (sampleLohEvents != null)
            mSampleLOHData.addAll(sampleLohEvents.stream().filter(x -> !x.Skipped).collect(Collectors.toList()));

        // Handle each of the 3 applicable types: DEL, BIALLELIC and AMP
        for (final DriverCatalog driverGene : mDriverCatalog)
        {
            HmfTranscriptRegion region = findRegion(driverGene.gene());

            if (region == null)
            {
                LOGGER.warn("driver gene({}) not found in all-genes map", driverGene.gene());
                continue;
            }

            final GeneCopyNumber geneCN = findGeneCopyNumber(driverGene);

            if(geneCN == null)
            {
                LOGGER.warn("gene({}) copy number data not found", driverGene.gene());
                continue;
            }

            DriverGeneData driverGeneData = new DriverGeneData(driverGene, region, geneCN);
            mDriverGeneDataList.add(driverGeneData);

            final List<SvBreakend> breakendList = mChrBreakendMap.get(region.chromosome());

            if (breakendList == null || breakendList.isEmpty())
            {
                writeDriverData(driverGeneData);
                continue;
            }

            if (driverGene.driver() == DriverType.DEL)
            {
                annotateDeleteEvent(driverGeneData, breakendList);
            }
            else if (driverGene.driver() == DriverType.AMP)
            {
                annotateAmplification(driverGeneData, breakendList);
            }
            else
            {
                // treat DNDS and HOTSPOT as potentially biallelic events
                annotateBiallelicEvent(driverGeneData);
            }
        }

        mSampleLOHData.clear();
        mChrBreakendMap = null;

        mPerfCounter.stop();
    }

    private void annotateDeleteEvent(final DriverGeneData driverGeneData, final List<SvBreakend> breakendList)
    {
        /* DEL identification:
            - 1 or 2 SVs which caused this, start from DEL region (ie gene) and work out
            - gene copy number table - minCopyRegion < 0.5 makes it a DEL-type driver candidate
            - walk out from here to end of LOH event ie where (1- BAF*CN) > 0.5 again
            - see COLO829T and PTEN as an example - here there is a simple DEL and then most of the other chromatid has been lost
         */

        // first find the min copy number within the gene region
        // then walk out in both directions to find the SV which caused the loss
        // and then walk out again until heterozygosity is gained (ie the end of an LOH)

        SvBreakend minBreakend = null; // breakend with lowest copy number covering any part of the gene region
        boolean isStartBreakend = true;

        final DriverCatalog driverGene = driverGeneData.DriverGene;
        final GeneCopyNumber geneCN = driverGeneData.GeneCN;

        for (int i = 0; i < breakendList.size(); ++i)
        {
            final SvBreakend breakend = breakendList.get(i);
            final SvBreakend nextBreakend = i < breakendList.size() - 1 ? breakendList.get(i + 1) : null;

            // find the 2 breakends straddling the start of the gene region, and compare their CNs in the -1 / upwards direction
            if ((breakend.position() < geneCN.minRegionStart() && nextBreakend != null && nextBreakend.position() > geneCN.minRegionStart())
            || (minBreakend == null && breakend.position() > geneCN.minRegionStart() && breakend.position() < geneCN.minRegionEnd()))
            {
                // take the lower of the copy numbers (upwards-facing)
                if (nextBreakend == null || breakend.getCopyNumber(false) < nextBreakend.getCopyNumber(false)
                || nextBreakend.position() >= geneCN.minRegionEnd())
                {
                    minBreakend = breakend;
                }
                else
                {
                    minBreakend = nextBreakend;
                }
            }
            else if(minBreakend == null && nextBreakend == null)
            {
                // no more breakends to consider - does this breakend definitely contribute to the DEL?
                minBreakend = breakend;
            }
            else if (breakend.position() >= geneCN.minRegionEnd())
            {
                if(minBreakend != null)
                    break;

                // some event beyond the gene caused loss so need to continue on until it's found

                // for look for an LOH event starting from here
                for (int j = i; j < breakendList.size(); ++j)
                {
                    final SvBreakend lohBreakend = breakendList.get(j);
                    if (isLOHEvent(lohBreakend, false))
                    {
                        minBreakend = lohBreakend;
                        isStartBreakend = false;
                        break;
                    }
                }

                if(minBreakend != null)
                    break;

                // otherwise select the first breakend facing away
                if(breakend.orientation() == -1)
                {
                    minBreakend = breakend;
                    break;
                }
            }
            else if (breakend.position() > geneCN.minRegionStart())
            {
                // check for a lower CN from another breakend inside the gene region
                if (breakend.getCopyNumber(false) < minBreakend.getCopyNumber(false))
                    minBreakend = breakend;
            }
        }

        if(minBreakend == null)
        {
            LOGGER.debug("sample({}) gene({}) not allocated to SVs", mSampleId, geneToStr(driverGene, driverGeneData.Region));
            writeDriverData(driverGeneData);
            return;
        }

        LOGGER.debug(String.format("gene(%s) min copy number at breakend(%s cn=%.2f cnChg=%.2f)",
                driverGene.gene(), minBreakend.toString(), minBreakend.getCopyNumber(false),
                minBreakend.getSV().copyNumberChange(minBreakend.usesStart())));

        SvBreakend startBreakend = isStartBreakend ? minBreakend : null;
        SvBreakend endBreakend = !isStartBreakend ? minBreakend : null;

        if(startBreakend != null)
        {
            // now walk forwards and backwards to find the next SVs
            int startIndex = startBreakend.getChrPosIndex();
            endBreakend = findDeletionBreakend(breakendList, startIndex, true, false);

            driverGeneData.addSvBreakend(startBreakend, "MIN");
        }

        if(endBreakend != null)
            driverGeneData.addSvBreakend(endBreakend, "MIN");

        // find straddling LOH for this DEL
        SvBreakend preStartBreakend = null;
        SvBreakend postEndBreakend = null;

        SvLOH matchedLohEvent = null;
        if(startBreakend != null && endBreakend != null)
        {
            for (final SvLOH lohEvent : mSampleLOHData)
            {
                if(lohEvent.Chromosome.equals(startBreakend.chromosome())
                && lohEvent.PosStart <= startBreakend.position() && lohEvent.PosEnd >= endBreakend.position())
                {
                    matchedLohEvent = lohEvent;
                    preStartBreakend = lohEvent.getBreakend(true);
                    postEndBreakend = lohEvent.getBreakend(false);
                    break;
                }
            }
        }

        // look to next event
        if(preStartBreakend != null && postEndBreakend != null)
        {
            driverGeneData.addSvBreakend(preStartBreakend, "LOH");
            driverGeneData.addSvBreakend(postEndBreakend, "LOH");
        }
        else if(preStartBreakend != null)
        {
            driverGeneData.addSvBreakend(preStartBreakend, "LOH_" + matchedLohEvent.SegEnd);
        }
        else if(postEndBreakend != null)
        {
            driverGeneData.addSvBreakend(postEndBreakend, "LOH_" + matchedLohEvent.SegStart);
        }
        else
        {
            driverGeneData.setMissedLohSVs(true);
        }

        writeDriverData(driverGeneData);
    }

    private SvBreakend findDeletionBreakend(final List<SvBreakend> breakendList, int startIndex, boolean walkForwards, boolean requireGOH)
    {
        int index = walkForwards ? startIndex + 1 : startIndex - 1;

        while(index >= 0 && index <= breakendList.size() - 1)
        {
            final SvBreakend breakend = breakendList.get(index);

            // this first variant following the DEL region ought to be facing away
            if ((walkForwards && breakend.orientation() == -1) || (!walkForwards && breakend.orientation() == 1))
            {
                if(requireGOH)
                {
                    if(isLOHEvent(breakend, !walkForwards))
                        return breakend;
                }
                else
                {
                    return breakend;
                }
            }

            index = walkForwards ? index + 1 : index - 1;
        }

        return null;
    }

    private boolean isLOHEvent(final SvBreakend breakend, boolean checkStart)
    {
        for (final SvLOH lohEvent : mSampleLOHData)
        {
            if ((checkStart && lohEvent.getBreakend(true) == breakend)
            || (!checkStart && lohEvent.getBreakend(false) == breakend))
            {
                return true;
            }
        }

        return false;
    }

    private void annotateBiallelicEvent(DriverGeneData driverGeneData)
    {
        final DriverCatalog driverGene = driverGeneData.DriverGene;
        final GeneCopyNumber geneCN = driverGeneData.GeneCN;

        // for biallelic events, find the straddling LOH event
        for (final SvLOH lohEvent : mSampleLOHData)
        {
            if(!lohEvent.Chromosome.equals(geneCN.chromosome()))
                continue;

            if(lohEvent.PosStart > geneCN.minRegionEnd() || lohEvent.PosEnd < geneCN.minRegionStart())
                continue;

            // now find the corresponding breakends
            SvBreakend startBreakend = lohEvent.getBreakend(true);
            SvBreakend endBreakend = lohEvent.getBreakend(false);

            if(startBreakend == null && endBreakend == null)
            {
                LOGGER.debug("sample({}) gene({}) LOH not found for driver gene",
                        mSampleId, geneToStr(driverGene, driverGeneData.Region));
                break;
            }


            if(startBreakend != null && endBreakend != null)
            {
                driverGeneData.addSvBreakend(startBreakend, "LOH");
                driverGeneData.addSvBreakend(endBreakend, "LOH");
            }
            else if(startBreakend != null)
            {
                driverGeneData.addSvBreakend(startBreakend, "LOH_" + lohEvent.SegEnd);
            }
            else if(endBreakend != null)
            {
                driverGeneData.addSvBreakend(endBreakend, "LOH_" + lohEvent.SegStart);
            }

            writeDriverData(driverGeneData);
            return;
        }

        driverGeneData.setMissedLohSVs(true);
        writeDriverData(driverGeneData);
    }

    private void annotateAmplification(final DriverGeneData driverGeneData, final List<SvBreakend> breakendList)
    {
        // find the cause - DUP, foldback, otherwise assume whole-chromatid duplication

        // trying to find the breakends which amplify this gene
        // take any INV, DUP or TI which straddles the gene

        final DriverCatalog driverGene = driverGeneData.DriverGene;
        final HmfTranscriptRegion region = driverGeneData.Region;

        double maxCopyNumber = 0;
        SvVarData maxSvStart = null;
        SvVarData maxSvEnd = null;
        SvBreakend foldbackBreakend = null;

        for (int i = 0; i < breakendList.size(); ++i)
        {
            final SvBreakend breakend = breakendList.get(i);

            if(breakend.position() > region.start())
            {
                // past the gene but AMP may be explained by an unchained foldback
                if(foldbackBreakend == null && breakend.orientation() == 1
                && !breakend.getSV().getFoldbackLink(breakend.usesStart()).isEmpty())
                {
                    foldbackBreakend = breakend;
                    break;
                }

                continue; // looking for a foldback as an explanation
            }

            final SvVarData varStart = breakend.getSV();

            if(varStart.getCluster().getResolvedType() == RESOLVED_TYPE_LOW_QUALITY)
                continue;

            if(breakend.orientation() == -1 && !varStart.getFoldbackLink(breakend.usesStart()).isEmpty())
                foldbackBreakend = breakend;

            if(varStart.type() == DUP || varStart.type() == INV)
            {
                if(varStart.position(true) <= region.start() && varStart.position(false) >= region.end())
                {
                    LOGGER.debug(String.format("sample(%s) cluster(%s) gene(%s) single SV(%s %s) cn(%.2f) cnChg(%.2f)",
                            mSampleId, varStart.getCluster().id(), geneToStr(driverGene, region), varStart.posId(), varStart.type(),
                            varStart.copyNumber(true), varStart.copyNumberChange(true)));

                    driverGeneData.addSvBreakend(varStart.getBreakend(true), "SV");
                    driverGeneData.addSvBreakend(varStart.getBreakend(false), "SV");

                    if(varStart.copyNumberChange(true) > maxCopyNumber)
                    {
                        maxCopyNumber = varStart.copyNumberChange(true);
                        maxSvStart = varStart;
                        maxSvEnd = varStart;
                    }

                    continue;
                }
            }

            final SvLinkedPair tiPair = varStart.getLinkedPair(breakend.usesStart());

            if(tiPair == null)
                continue;

            if(breakend.position() + tiPair.length() < region.end())
                continue;

            final SvVarData varEnd = tiPair.getOtherSV(varStart);
            final SvBreakend beStart = tiPair.getBreakend(true);
            final SvBreakend beEnd = tiPair.getBreakend(false);

            LOGGER.debug(String.format("sample(%s) cluster(%d fb=%s) gene(%s) SVs start(%s cn=%.2f cnChg=%.2f) end(%s cn=%.2f cnChg=%.2f) in linked pair",
                    mSampleId, varStart.getCluster().id(), varStart.getCluster().getFoldbacks().size(), geneToStr(driverGene, region),
                    varStart.posId(), varStart.copyNumber(beStart.usesStart()), varStart.copyNumberChange(beStart.usesStart()),
                    varEnd.posId(), varEnd.copyNumber(beEnd.usesStart()), varEnd.copyNumberChange(beEnd.usesStart())));

            driverGeneData.addSvBreakend(beStart, "TI");
            driverGeneData.addSvBreakend(beEnd, "TI");

            if(varStart.copyNumberChange(true) > maxCopyNumber)
            {
                maxCopyNumber = varStart.copyNumberChange(true);
                maxSvStart = varStart;
                maxSvEnd = varEnd;
            }
        }

        if(maxSvStart == null || maxSvEnd == null)
        {
            if(foldbackBreakend != null)
                driverGeneData.addSvBreakend(foldbackBreakend, "FB");
        }

        writeDriverData(driverGeneData);
    }

    private static final String geneToStr(DriverCatalog driverGene, HmfTranscriptRegion region)
    {
        return String.format("%s: %s %s:%d-%d",
                driverGene.driver(), driverGene.gene(), region.chromosome(), region.start(), region.end());
    }

    private final SvLOH findLOHEventForRegion(final GeneCopyNumber geneCN)
    {
        // check for any LOH not linked to an SV which straddles this region
        for (final SvLOH lohEvent : mSampleLOHData)
        {
            if (lohEvent.Chromosome.equals(geneCN.chromosome())
            && lohEvent.PosStart <= geneCN.minRegionStart() && lohEvent.PosEnd >= geneCN.minRegionEnd())
            {
                return lohEvent;
            }
        }

        return null;
    }

    private final GeneCopyNumber findGeneCopyNumber(final DriverCatalog driverGene)
    {
        for(GeneCopyNumber geneCN : mGeneCopyNumberData)
        {
            if(driverGene.gene().equals(geneCN.gene()))
                return geneCN;
        }

        return null;
    }

    private void writeDriverData(final DriverGeneData driverGeneData)
    {
        try
        {
            if(mFileWriter == null)
            {
                String outputFileName = mOutputDir;

                outputFileName += "SVA_DRIVERS.csv";

                mFileWriter = createBufferedWriter(outputFileName, false);

                mFileWriter.write("SampleId,Gene,GeneType,DriverType,DriverLikelihood");
                mFileWriter.write(",ClusterId,SvId,IsStart,MatchInfo");
                mFileWriter.write(",SamplePloidy,Chromosome,Arm,MinCN,CentromereCN,TelomereCN");
                mFileWriter.newLine();
            }

            BufferedWriter writer = mFileWriter;

            final DriverCatalog driverGene = driverGeneData.DriverGene;
            final HmfTranscriptRegion region = driverGeneData.Region;
            double[] cnData = mChrCopyNumberMap.get(region.chromosome());
            final String arm = getChromosomalArm(region.chromosome(), region.geneStart());
            final GeneCopyNumber geneCN = driverGeneData.GeneCN;
            int refClusterId = -1;

            final List<SvBreakend> svBreakends = driverGeneData.getSvBreakends();

            int dataCount;

            if(svBreakends.isEmpty())
            {
                dataCount = 1;
            }
            else
            {
                dataCount = svBreakends.size();

                if(driverGeneData.missedLohSVs())
                    ++dataCount;
            }

            for(int i = 0; i < dataCount; ++i)
            {
                writer.write(String.format("%s,%s,%s,%s,%.4f",
                        mSampleId, driverGene.gene(), driverGene.category(), driverGene.driver(), driverGene.driverLikelihood()));

                if(svBreakends.isEmpty() || i >= svBreakends.size())
                {
                    final SvLOH lohEvent = driverGene.driver() != DriverType.AMP ? findLOHEventForRegion(geneCN) : null;

                    if (lohEvent != null)
                        writer.write(String.format(",-1,,,%s;%s", lohEvent.SegStart, lohEvent.SegEnd));
                    else
                        writer.write(String.format(",-1,,,"));
                }
                else
                {
                    final SvBreakend breakend = svBreakends.get(i);
                    final String svInfo = driverGeneData.getSvInfoList().get(i);

                    if(breakend == null)
                        break;

                    final SvVarData var = breakend.getSV();
                    refClusterId = var.getCluster().id();

                    // cache info against the SV
                    var.setDriveGene(String.format("%s;%s;%s",
                            driverGeneData.DriverGene.driver(), driverGeneData.DriverGene.gene(), svInfo), breakend.usesStart());

                    writer.write(String.format(",%d,%s,%s,%s", var.getCluster().id(), var.origId(), breakend.usesStart(), svInfo));
                }

                writer.write(String.format(",%.2f,%s,%s,%.2f,%.2f,%.2f",
                        mSamplePloidy, region.chromosome(), arm, geneCN != null ? geneCN.minCopyNumber() : -1,
                        cnData[CENTROMERE_CN],
                        arm == CHROMOSOME_ARM_P ? cnData[P_ARM_TELOMERE_CN] : cnData[Q_ARM_TELOMERE_CN]));

                writer.newLine();
            }

            mVisWriter.addGeneExonData(refClusterId, region.geneID(), region.gene(),
                    "", region.chromosome(), "DRIVER");

        }
        catch (final IOException e)
        {
            LOGGER.error("error writing cluster-data to outputFile: {}", e.toString());
        }
    }

    public void close()
    {
        mPerfCounter.logStats();

        closeBufferedWriter(mFileWriter);
        closeBufferedWriter(mGCNFileWriter);
    }

    private void checkPseudoGeneAnnotations(final List<SvCluster> clusters)
    {
        for(final SvCluster cluster : clusters)
        {
            // isSpecificCluster(cluster);

            GeneAnnotation pseudoGene = null;
            String transcriptId = "";

            for(final SvLinkedPair pair : cluster.getLinkedPairs())
            {
                if(pair.length() > SHORT_TI_LENGTH * 8)
                    continue;

                final SvBreakend lower = pair.getBreakend(true);
                final SvBreakend upper = pair.getBreakend(false);

                // for any TI falling within the same gene, check for an exon boundary match
                if(lower.getSV().getGenesList(lower.usesStart()).isEmpty() || upper.getSV().getGenesList(upper.usesStart()).isEmpty())
                    continue;

                for(final GeneAnnotation gene1 : lower.getSV().getGenesList(lower.usesStart()))
                {
                    for(final GeneAnnotation gene2 : upper.getSV().getGenesList(upper.usesStart()))
                    {
                        if(!gene1.GeneName.equals(gene2.GeneName))
                            continue;

                        final String exonData[] = mGeneTranscriptCollection.getExonDetailsForPosition(gene1, lower.position(), upper.position());

                        if(exonData[PSEUDO_GENE_DATA_TRANS_ID] != null)
                        {
                            pseudoGene = gene1;
                            transcriptId = exonData[PSEUDO_GENE_DATA_TRANS_ID];

                            String exonMatchData = String.format("%s;%s;%s;%s",
                                    transcriptId, exonData[PSEUDO_GENE_DATA_EXON_RANK],
                                    exonData[PSEUDO_GENE_DATA_EXON_MAX], exonData[PSEUDO_GENE_DATA_EXON_LENGTH]);


                            pair.setExonMatchData(exonMatchData);
                        }
                    }
                }
            }

            if(pseudoGene != null)
            {
                mVisWriter.addGeneExonData(cluster.id(), pseudoGene.StableId, pseudoGene.GeneName,
                        transcriptId, pseudoGene.chromosome(), "PSEUDO");
            }
        }
    }

    private static String GENE_CN_DATA_FILE = "SVA_GENE_COPY_NUMBER.csv";
    private static int GENE_CN_DATA_FILE_ITEM_COUNT = 30;

    private void loadGeneCopyNumberDataFile(final String gcnFileName)
    {
        if(gcnFileName.isEmpty() || !Files.exists(Paths.get(gcnFileName)))
            return;

        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(gcnFileName));

            // skip field names
            String line = fileReader.readLine();

            if (line == null)
            {
                LOGGER.error("Empty gene copy number CSV file({})", gcnFileName);
                return;
            }

            String currentSample = "";
            List<GeneCopyNumber> gcnDataList = null;

            while ((line = fileReader.readLine()) != null)
            {
                String[] items = line.split(",");

                if (items.length != GENE_CN_DATA_FILE_ITEM_COUNT)
                    continue;

                String sampleId = items[1];

                if(currentSample.isEmpty() || !currentSample.equals(sampleId))
                {
                    gcnDataList = Lists.newArrayList();
                    currentSample = sampleId;
                    mSampleGeneCopyNumberMap.put(currentSample, gcnDataList);

                }

                //                          0       1       2           3   4   5      6              7           8
//                mGCNFileWriter.write("modified,sampleId,chromosome,start,end,gene,chromosomeBand,transcriptId,transcriptVersion");
                //                       10             11           12             13                  14                 15        16              17             18                  19                      20
//                mGCNFileWriter.write(",minCopyNumber,maxCopyNumber,somaticRegions,germlineHomRegions,germlineHetRegions,minRegions,minRegionStart,minRegionEnd,minRegionStartSupport,minRegionEndSupport,minRegionMethod");
                //                          21
//                mGCNFileWriter.write(",nonsenseBiallelicVariants,nonsenseNonBiallelicVariants,nonsenseNonBiallelicPloidy,spliceBiallelicVariants,spliceNonBiallelicVariants,spliceNonBiallelicPloidy");
//                mGCNFileWriter.write(",missenseBiallelicVariants,missenseNonBiallelicVariants,missenseNonBiallelicPloidy,minMinorAllelePloidy");

                int index = 2;

                GeneCopyNumber gcnData = ImmutableGeneCopyNumber.builder()
                        .chromosome(items[index++])
                        .start(Integer.parseInt(items[index++]))
                        .end(Integer.parseInt(items[index++]))
                        .gene(items[index++])
                        .chromosomeBand(items[index++])
                        .transcriptID(items[index++])
                        .transcriptVersion(Integer.parseInt(items[index++]))
                        .minCopyNumber(Double.parseDouble(items[index++]))
                        .maxCopyNumber(Double.parseDouble(items[index++]))
                        .somaticRegions(Integer.parseInt(items[index++]))
                        .germlineHomRegions(Integer.parseInt(items[index++]))
                        .germlineHet2HomRegions(Integer.parseInt(items[index++]))
                        .minRegions(Integer.parseInt(items[index++]))
                        .minRegionStart(Integer.parseInt(items[index++]))
                        .minRegionEnd(Integer.parseInt(items[index++]))
                        .minRegionStartSupport(SegmentSupport.valueOf(items[index++]))
                        .minRegionEndSupport(SegmentSupport.valueOf(items[index++]))
                        .minRegionMethod(CopyNumberMethod.valueOf(items[index++]))
                        .nonsenseBiallelicCount(Integer.parseInt(items[index++]))
                        .nonsenseNonBiallelicCount(Integer.parseInt(items[index++]))
                        .nonsenseNonBiallelicPloidy(Double.parseDouble(items[index++]))
                        .spliceBiallelicCount(Integer.parseInt(items[index++]))
                        .spliceNonBiallelicCount(Integer.parseInt(items[index++]))
                        .spliceNonBiallelicPloidy(Double.parseDouble(items[index++]))
                        .missenseBiallelicCount(Integer.parseInt(items[index++]))
                        .missenseNonBiallelicCount(Integer.parseInt(items[index++]))
                        .missenseNonBiallelicPloidy(Double.parseDouble(items[index++]))
                        .minMinorAllelePloidy(Double.parseDouble(items[index++]))
                        .build();


                gcnDataList.add(gcnData);
            }
        }
        catch (IOException e)
        {
            LOGGER.error("Failed to read gene copy number CSV file({}): {}", gcnFileName, e.toString());
        }
    }

    private void writeGeneCopyNumberData(final GeneCopyNumber gcnData)
    {
        try
        {
            if(mGCNFileWriter == null)
            {
                String outputFileName = mOutputDir;

                outputFileName += GENE_CN_DATA_FILE;

                mGCNFileWriter = createBufferedWriter(outputFileName, false);

                mGCNFileWriter.write("modified,sampleId,chromosome,start,end,gene,chromosomeBand,transcriptId,transcriptVersion");
                mGCNFileWriter.write(",minCopyNumber,maxCopyNumber,somaticRegions,germlineHomRegions,germlineHetRegions,minRegions,minRegionStart,minRegionEnd,minRegionStartSupport,minRegionEndSupport,minRegionMethod");
                mGCNFileWriter.write(",nonsenseBiallelicVariants,nonsenseNonBiallelicVariants,nonsenseNonBiallelicPloidy,spliceBiallelicVariants,spliceNonBiallelicVariants,spliceNonBiallelicPloidy");
                mGCNFileWriter.write(",missenseBiallelicVariants,missenseNonBiallelicVariants,missenseNonBiallelicPloidy,minMinorAllelePloidy");
                mGCNFileWriter.newLine();
            }

            BufferedWriter writer = mGCNFileWriter;

            final Timestamp timestamp = new Timestamp(new Date().getTime());

            writer.write(String.format("%s,%s,%s,%d,%d,%s,%s,%s,%d",
                    timestamp, mSampleId, gcnData.chromosome(), gcnData.start(), gcnData.end(),
                    gcnData.gene(), gcnData.chromosomeBand(), gcnData.transcriptID(), gcnData.transcriptVersion()));

            writer.write(String.format(",%.4f,%.4f,%d,%d,%d,%d,%d,%d,%s,%s,%s",
                    gcnData.minCopyNumber(), gcnData.maxCopyNumber(), gcnData.somaticRegions(),
                    gcnData.germlineHomRegions(), gcnData.germlineHet2HomRegions(),
                    gcnData.minRegions(), gcnData.minRegionStart(), gcnData.minRegionEnd(),
                    gcnData.minRegionStartSupport(), gcnData.minRegionEndSupport(), gcnData.minRegionMethod().toString()));


            writer.write(String.format(",%d,%d,%.4f,%d,%d,%.4f,%d,%d,%.4f,%.4f,",
                    gcnData.nonsenseBiallelicCount(), gcnData.nonsenseNonBiallelicCount(), gcnData.nonsenseNonBiallelicPloidy(),
                    gcnData.spliceBiallelicCount(), gcnData.spliceNonBiallelicCount(), gcnData.spliceNonBiallelicPloidy(),
                    gcnData.missenseBiallelicCount(), gcnData.missenseNonBiallelicCount(), gcnData.missenseNonBiallelicPloidy(),
                    gcnData.minMinorAllelePloidy()));

            writer.newLine();

        }
        catch (final IOException e)
        {
            LOGGER.error("error writing gene copy number to outputFile: {}", e.toString());
        }
    }
}
