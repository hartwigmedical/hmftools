package com.hartwig.hmftools.linx.visualiser.file;

import static com.hartwig.hmftools.common.immune.ImmuneRegions.getIgRegion;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_SAMPLE_ID;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.linx.LinxConfig.LNX_LOGGER;
import static com.hartwig.hmftools.linx.LinxConfig.REF_GENOME_VERSION;
import static com.hartwig.hmftools.linx.analysis.ClusterClassification.isFilteredResolvedType;
import static com.hartwig.hmftools.common.purple.ChromosomeArm.P_ARM;
import static com.hartwig.hmftools.common.purple.ChromosomeArm.Q_ARM;
import static com.hartwig.hmftools.linx.types.SvBreakend.DIRECTION_CENTROMERE;
import static com.hartwig.hmftools.linx.types.SvBreakend.DIRECTION_TELOMERE;
import static com.hartwig.hmftools.linx.visualiser.file.VisGeneAnnotationType.EXON_LOST;
import static com.hartwig.hmftools.linx.visualiser.file.VisGeneAnnotationType.FUSION;
import static com.hartwig.hmftools.linx.visualiser.file.VisProteinDomain.PD_FIVE_PRIME_UTR;
import static com.hartwig.hmftools.linx.visualiser.file.VisProteinDomain.PD_NON_CODING;
import static com.hartwig.hmftools.linx.visualiser.file.VisProteinDomain.PD_THREE_PRIME_UTR;
import static com.hartwig.hmftools.linx.visualiser.file.VisSvData.INFO_TYPE_FOLDBACK;
import static com.hartwig.hmftools.linx.visualiser.file.VisSvData.INFO_TYPE_NORMAL;

import java.io.BufferedWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.gene.TranscriptProteinData;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.linx.chaining.SvChain;
import com.hartwig.hmftools.linx.cn.SvCNData;
import com.hartwig.hmftools.linx.types.LinkedPair;
import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.types.SvVarData;

public class VisDataWriter
{
    // references only
    private EnsemblDataCache mGeneDataCache;

    private boolean mEnabled;
    private final String mOutputDir;
    private final boolean mGermline; // currently incomplete and disabled

    private boolean mBatchOutput;
    private BufferedWriter mSvFileWriter;
    private BufferedWriter mSegmentFileWriter;
    private BufferedWriter mCnFileWriter;
    private BufferedWriter mGeneFileWriter;
    private BufferedWriter mProteinDomainFileWriter;
    private BufferedWriter mFusionFileWriter;

    public static final String COHORT_VIS_SVS_FILE = "LNX_VIS_SVS.tsv";
    public static final String COHORT_VIS_LINKS_FILE = "LNX_VIS_SEGMENTS.tsv";
    public static final String COHORT_VIS_COPY_NUMBER_FILE = "LNX_VIS_COPY_NUMBER.tsv";
    public static final String COHORT_VIS_GENE_EXONS_FILE = "LNX_VIS_GENE_EXONS.tsv";
    public static final String COHORT_VIS_PROTEIN_FILE = "LNX_VIS_PROTEIN_DOMAINS.tsv";
    public static final String COHORT_VIS_FUSIONS_FILE = "LNX_VIS_FUSIONS.tsv";

    public static final DecimalFormat FORMAT = new DecimalFormat("0.0000");

    public VisDataWriter(
            final String outputDir, final EnsemblDataCache geneDataCache, boolean enabled, boolean isBatchOutput, boolean isGermline)
    {
        mEnabled = enabled;
        mGermline = isGermline;
        mOutputDir = outputDir;
        mBatchOutput = isBatchOutput;

        if(mBatchOutput && mEnabled)
        {
            initialiseCohortOutputFiles();
        }

        mGeneDataCache = geneDataCache;
    }

    private void initialiseCohortOutputFiles()
    {
        try
        {
            String sampleIdColumn = FLD_SAMPLE_ID + TSV_DELIM;
            mSvFileWriter = createBufferedWriter(mOutputDir + COHORT_VIS_SVS_FILE, false);
            mSvFileWriter.write(sampleIdColumn);
            mSvFileWriter.write(VisSvData.header());
            mSvFileWriter.newLine();

            mSegmentFileWriter = createBufferedWriter(mOutputDir + COHORT_VIS_LINKS_FILE, false);
            mSegmentFileWriter.write(sampleIdColumn);
            mSegmentFileWriter.write(VisSegment.header());
            mSegmentFileWriter.newLine();

            mCnFileWriter = createBufferedWriter(mOutputDir + COHORT_VIS_COPY_NUMBER_FILE, false);
            mCnFileWriter.write(sampleIdColumn);
            mCnFileWriter.write(VisCopyNumber.header());
            mCnFileWriter.newLine();

            mGeneFileWriter = createBufferedWriter(mOutputDir + COHORT_VIS_GENE_EXONS_FILE, false);
            mGeneFileWriter.write(sampleIdColumn);
            mGeneFileWriter.write(VisGeneExon.header());
            mGeneFileWriter.newLine();

            mProteinDomainFileWriter = createBufferedWriter(mOutputDir + COHORT_VIS_PROTEIN_FILE, false);
            mProteinDomainFileWriter.write(sampleIdColumn);
            mProteinDomainFileWriter.write(VisProteinDomain.header());
            mProteinDomainFileWriter.newLine();

            mFusionFileWriter = createBufferedWriter(mOutputDir + COHORT_VIS_FUSIONS_FILE, false);
            mFusionFileWriter.write(sampleIdColumn);
            mFusionFileWriter.write(VisFusion.header());
            mFusionFileWriter.newLine();

        }
        catch(IOException e)
        {
            LNX_LOGGER.error("failed to open and write output file headers");
        }
    }

    public void close()
    {
        closeBufferedWriter(mSvFileWriter);
        closeBufferedWriter(mSegmentFileWriter);
        closeBufferedWriter(mCnFileWriter);
        closeBufferedWriter(mGeneFileWriter);
        closeBufferedWriter(mProteinDomainFileWriter);
        closeBufferedWriter(mFusionFileWriter);
    }

    public synchronized void writeOutput(
            final VisSampleData sampleData, final List<SvCluster> clusters, final List<SvVarData> variants,
            final Map<String,List<SvCNData>> chrCnDataMap)
    {
        if(!mEnabled)
            return;

        writeSvData(sampleData, variants);
        writeSegmentData(sampleData, clusters);
        writeGeneData(sampleData);
        writeCopyNumberData(sampleData, chrCnDataMap);
        writeFusions(sampleData);
    }

    private void writeSvData(final VisSampleData sampleData, final List<SvVarData> variants)
    {
        List<VisSvData> svDataList = Lists.newArrayList();

        for(final SvVarData var : variants)
        {
            final SvCluster cluster = var.getCluster();

            final List<SvChain> chains = cluster.findChains(var);

            // repeat an SV for every time it appears in a chain
            int chainCount = chains.isEmpty() ? 1 : chains.size();
            int unchainedChainId = chains.isEmpty() ? cluster.getChainId(var) : -1;

            for(int i = 0; i < chainCount; ++i)
            {
                final SvChain chain = !chains.isEmpty() ? chains.get(i) : null;
                int chainId = chain == null ? unchainedChainId : chain.id();

                final SvBreakend beStart = var.getBreakend(true);
                final SvBreakend beEnd = var.getBreakend(false);

                svDataList.add(new VisSvData(sampleData.sampleId(), cluster.id(), chainId, var.id(),
                        var.type(), cluster.getResolvedType(), cluster.isSyntheticType(),
                        beStart.chromosome(), beEnd != null ? beEnd.chromosome() : "-1",
                        beStart.position(), beEnd != null ? beEnd.position() : 0,
                        beStart.orientation(), beEnd != null ? beEnd.orientation() : 0,
                        beStart.getSV().getFoldbackBreakend(beStart.usesStart()) == null ? INFO_TYPE_NORMAL : INFO_TYPE_FOLDBACK,
                        beEnd!= null ? (beEnd.getSV().getFoldbackBreakend(beEnd.usesStart()) == null ? INFO_TYPE_NORMAL : INFO_TYPE_FOLDBACK) : "",
                        var.jcn(), chain != null && chain.isDoubleMinute()));
            }
        }

        try
        {
            if(mBatchOutput)
            {
                for(final VisSvData data : svDataList)
                {
                    mSvFileWriter.write(sampleData.sampleId() + TSV_DELIM);
                    mSvFileWriter.write(VisSvData.toString(data));
                    mSvFileWriter.newLine();
                }
            }
            else
            {
                VisSvData.write(VisSvData.generateFilename(mOutputDir, sampleData.sampleId(), mGermline), svDataList);
            }
        }
        catch(IOException e)
        {
            LNX_LOGGER.error("filed to write VIS SV data output");
        }
    }

    private void writeSegmentData(final VisSampleData sampleData, final List<SvCluster> clusters)
    {
        // write out the links from each chain and a link from the chain-end breakends to telomere or centromere
        List<VisSegment> segments = Lists.newArrayList();

        for(final SvCluster cluster : clusters)
        {
            if(cluster.getSvCount() == 1)
            {
                if(isFilteredResolvedType(cluster.getResolvedType()))
                    continue;

                // BNDs, INVs and SGLs will have lines showing their orientation
                // simple SVs will show connectors when plotting a single cluster
            }

            // for any linked pair which is repeated in a separate chain, skip writing it for subsequent chains
            List<LinkedPair> uniquePairs = Lists.newArrayList();

            // log chains in order of highest to lowest ploidy so that where an SV is is more than 1 chain it will show the max ploidy
            List<SvChain> chains = Lists.newArrayList();

            for(final SvChain chain : cluster.getChains())
            {
                int index = 0;
                while(index < chains.size())
                {
                    if(chain.jcn() > chains.get(index).jcn())
                        break;
                    else
                        ++index;
                }

                chains.add(index, chain);
            }

            for(final SvChain chain : chains)
            {
                // log the start of the chain
                boolean startsOnEnd = chain.getFirstSV() == chain.getLastSV(); // closed loop chains eg DMs
                double chainPloidy = chain.jcn();

                if(!chain.isClosedLoop())
                {
                    SvBreakend breakend = chain.getOpenBreakend(true);

                    if(breakend != null)
                    {
                        segments.add(new VisSegment(sampleData.sampleId(), cluster.id(), chain.id(), breakend.chromosome(),
                                getPositionValue(breakend, true), getPositionValue(breakend, false), chainPloidy, false));
                    }
                }

                for(final LinkedPair pair : chain.getLinkedPairs())
                {
                    boolean isRepeat = uniquePairs.stream().anyMatch(x -> x.matches(pair));

                    if(isRepeat)
                        continue;

                    uniquePairs.add(pair);

                    double linkPloidy = 0;

                    if(cluster.requiresReplication() || cluster.getChains().size() > 1)
                    {
                        for(final SvChain otherChain : cluster.getChains())
                        {
                            int linkRepeats = (int) otherChain.getLinkedPairs().stream().filter(x -> x.matches(pair)).count();
                            linkPloidy += linkRepeats * otherChain.jcn();
                        }
                    }
                    else
                    {
                        linkPloidy = chainPloidy;
                    }

                    final SvBreakend beStart = pair.getBreakend(true);
                    final SvBreakend beEnd = pair.getBreakend(false);

                    segments.add(new VisSegment(sampleData.sampleId(), cluster.id(), chain.id(), beStart.chromosome(),
                            String.valueOf(beStart.position()), String.valueOf(beEnd.position()), linkPloidy, chain.isDoubleMinute()));
                }

                if(!chain.isClosedLoop())
                {
                    // log the end of the chain out to centromere or telomere
                    SvBreakend breakend = chain.getOpenBreakend(false);

                    if(breakend != null && !startsOnEnd)
                    {
                        segments.add(new VisSegment(sampleData.sampleId(), cluster.id(), chain.id(), breakend.chromosome(),
                                getPositionValue(breakend, true), getPositionValue(breakend, false), chainPloidy, false));
                    }
                }
            }

            // finally write out all unchained SVs with telomere and centromere links shown
            for(final SvVarData var : cluster.getUnlinkedSVs())
            {
                int chainId = cluster.getChainId(var);

                for(int be = SE_START; be <= SE_END; ++be)
                {
                    final SvBreakend breakend = var.getBreakend(be);

                    if(breakend == null)
                        continue;

                    segments.add(new VisSegment(sampleData.sampleId(), cluster.id(), chainId, breakend.chromosome(),
                            getPositionValue(breakend, true), getPositionValue(breakend, false), var.jcn(), false));
                }
            }
        }

        try
        {
            if(mBatchOutput)
            {
                for(final VisSegment data : segments)
                {
                    mSegmentFileWriter.write(sampleData.sampleId() + TSV_DELIM);
                    mSegmentFileWriter.write(VisSegment.toString(data));
                    mSegmentFileWriter.newLine();
                }
            }
            else
            {
                VisSegment.write(VisSegment.generateFilename(mOutputDir, sampleData.sampleId(), mGermline), segments);
            }
        }
        catch (final IOException e)
        {
            LNX_LOGGER.error("error writing to visualiser segments file: {}", e.toString());
        }
    }

    private String getPositionValue(final SvBreakend breakend, boolean isChainEnd)
    {
        int position = breakend.position();

        if(breakend.orientation() == 1 && breakend.arm() == P_ARM)
        {
            return isChainEnd ? DIRECTION_TELOMERE : Long.toString(position);
        }
        else if(breakend.orientation() == -1 && breakend.arm() == P_ARM)
        {
            return isChainEnd ? Long.toString(position) : DIRECTION_CENTROMERE;
        }
        else if(breakend.orientation() == -1 && breakend.arm() == Q_ARM)
        {
            return isChainEnd ? Long.toString(position) : DIRECTION_TELOMERE;
        }
        else
        {
            return isChainEnd ? DIRECTION_CENTROMERE : Long.toString(position);
        }
    }

    private void writeGeneData(final VisSampleData sampleData)
    {
        if(!mEnabled)
            return;

        final List<VisGeneExon> geneExonList = Lists.newArrayList();
        final List<VisProteinDomain> proteinList = Lists.newArrayList();

        for(final VisGeneData geneData : sampleData.getGeneData())
        {
            if(geneData.ClusterId < 0) // skip genes not linked to a cluster
                continue;

            if(checkAddIgExonRegions(sampleData.sampleId(), geneData, geneExonList))
                continue;

            final TranscriptData transData = mGeneDataCache.getTranscriptData(geneData.GeneId, geneData.TransName);

            if(transData == null || transData.exons().isEmpty())
                continue;

            for(final ExonData exonData : transData.exons())
            {
                if(geneData.ExonPositionOffsets.isEmpty())
                {
                    geneExonList.add(new VisGeneExon(sampleData.sampleId(), geneData.ClusterId, geneData.GeneName, transData.TransName,
                            geneData.Chromosome, geneData.AnnotationType, exonData.Rank, exonData.Start, exonData.End));
                }
                else
                {
                    // used to more accurately plot pseudogene gene deletions
                    final int[] exonPosOffsets = geneData.ExonPositionOffsets.get(exonData.Rank);
                    final int[] exonsLost = geneData.ExonsLostOffsets.get(exonData.Rank);

                    int exonStart = exonData.Start + exonPosOffsets[SE_START];
                    int exonEnd = exonData.End + exonPosOffsets[SE_END];

                    geneExonList.add(new VisGeneExon(sampleData.sampleId(), geneData.ClusterId, geneData.GeneName, transData.TransName,
                            geneData.Chromosome, geneData.AnnotationType, exonData.Rank, exonStart, exonEnd));

                    if(exonsLost != null)
                    {
                        geneExonList.add(new VisGeneExon(sampleData.sampleId(), geneData.ClusterId, geneData.GeneName, transData.TransName,
                                geneData.Chromosome, EXON_LOST, exonData.Rank,
                                exonStart + exonsLost[SE_START], exonEnd + exonsLost[SE_END]));
                    }
                }
            }

            int transId = geneData.TransId > 0 ? geneData.TransId : transData.TransId;

            final List<TranscriptProteinData> transProteinData = mGeneDataCache.getTranscriptProteinDataMap().get(transId);

            if(transProteinData != null)
            {
                for(final TranscriptProteinData proteinData : transProteinData)
                {
                    final Integer[] domainPositions = EnsemblDataCache.getProteinDomainPositions(proteinData, transData);

                    if(domainPositions[SE_START] != null && domainPositions[SE_END] != null)
                    {
                        proteinList.add(new VisProteinDomain(sampleData.sampleId(), geneData.ClusterId, transData.TransName, geneData.Chromosome,
                                domainPositions[SE_START], domainPositions[SE_END], proteinData.HitDescription));
                    }
                }
            }

            // show the 5' and 3' UTR or non-coding regions as 'protein domains'
            if(transData.CodingStart != null && transData.CodingEnd != null)
            {
                int fivePrimeUtrStart = transData.Strand == 1 ? transData.TransStart : transData.CodingEnd + 1;
                int fivePrimeUtrEnd = transData.Strand == 1 ? transData.CodingStart - 1 : transData.TransEnd;

                int threePrimeUtrStart = transData.Strand == 1 ? transData.CodingEnd + 1 : transData.TransStart;
                int threePrimeUtrEnd = transData.Strand == 1 ? transData.TransEnd : transData.CodingStart - 1;

                if(fivePrimeUtrStart < fivePrimeUtrEnd)
                {
                    proteinList.add(new VisProteinDomain(sampleData.sampleId(), geneData.ClusterId, transData.TransName, geneData.Chromosome,
                            fivePrimeUtrStart, fivePrimeUtrEnd, PD_FIVE_PRIME_UTR));
                }

                if(threePrimeUtrStart < threePrimeUtrEnd)
                {
                    proteinList.add(new VisProteinDomain(sampleData.sampleId(), geneData.ClusterId, transData.TransName, geneData.Chromosome,
                            threePrimeUtrStart, threePrimeUtrEnd, PD_THREE_PRIME_UTR));
                }
            }
            else
            {
                proteinList.add(new VisProteinDomain(sampleData.sampleId(), geneData.ClusterId, transData.TransName, geneData.Chromosome,
                        transData.TransStart, transData.TransEnd, PD_NON_CODING));
            }
        }

        try
        {
            if(mBatchOutput)
            {
                for(final VisGeneExon data : geneExonList)
                {
                    mGeneFileWriter.write(sampleData.sampleId() + TSV_DELIM);
                    mGeneFileWriter.write(VisGeneExon.toString(data));
                    mGeneFileWriter.newLine();
                }

                for(final VisProteinDomain data : proteinList)
                {
                    mProteinDomainFileWriter.write(sampleData.sampleId() + TSV_DELIM);
                    mProteinDomainFileWriter.write(VisProteinDomain.toString(data));
                    mProteinDomainFileWriter.newLine();
                }
            }
            else
            {
                VisGeneExon.write(VisGeneExon.generateFilename(mOutputDir, sampleData.sampleId(), mGermline), geneExonList);
                VisProteinDomain.write(VisProteinDomain.generateFilename(mOutputDir, sampleData.sampleId(), mGermline), proteinList);
            }
        }
        catch (final IOException e)
        {
            LNX_LOGGER.error("error writing to visualiser gene-exons file: {}", e.toString());
        }
    }

    private boolean checkAddIgExonRegions(final String sampleId, final VisGeneData geneData, final List<VisGeneExon> geneExonList)
    {
        if(!geneData.AnnotationType.equals(FUSION))
            return false;

        if(!geneData.TransName.contains("@IG"))
            return false;

        String igGene = geneData.TransName.replace("@", "");
        ChrBaseRegion igRegion = getIgRegion(igGene, REF_GENOME_VERSION);

        geneExonList.add(new VisGeneExon(sampleId, geneData.ClusterId, geneData.TransName, geneData.TransName,
                geneData.Chromosome, geneData.AnnotationType, 1, igRegion.start(), igRegion.end()));

        return true;
    }

    private void writeCopyNumberData(final VisSampleData sampleData, final Map<String,List<SvCNData>> chrCNDataMap)
    {
        List<VisCopyNumber> cnDataList = Lists.newArrayList();

        for(Map.Entry<String,List<SvCNData>> entry : chrCNDataMap.entrySet())
        {
            final String chromosome = entry.getKey();

            for(SvCNData cnData : entry.getValue())
            {
                cnDataList.add(new VisCopyNumber(
                        sampleData.sampleId(), chromosome, cnData.StartPos, cnData.EndPos, cnData.CopyNumber, cnData.ActualBaf));
            }
        }

        try
        {
            if(mBatchOutput)
            {
                for(final VisCopyNumber data : cnDataList)
                {
                    mCnFileWriter.write(sampleData.sampleId() + TSV_DELIM);
                    mCnFileWriter.write(VisCopyNumber.toString(data));
                    mCnFileWriter.newLine();
                }
            }
            else
            {
                VisCopyNumber.write(VisCopyNumber.generateFilename(mOutputDir, sampleData.sampleId(), mGermline), cnDataList);
            }
        }
        catch(IOException e)
        {
            LNX_LOGGER.error("filed to write VIS copy number output");
        }
    }

    private void writeFusions(final VisSampleData sampleData)
    {
        try
        {
            if(mBatchOutput)
            {
                for(final VisFusion visFusion : sampleData.getFusions())
                {
                    mFusionFileWriter.write(sampleData.sampleId() + TSV_DELIM);
                    mFusionFileWriter.write(VisFusion.toString(visFusion));
                    mFusionFileWriter.newLine();
                }
            }
            else
            {
                VisFusion.write(VisFusion.generateFilename(mOutputDir, sampleData.sampleId(), mGermline), sampleData.getFusions());
            }
        }
        catch(IOException e)
        {
            LNX_LOGGER.error("failed to write fusions vis file: {}", e.toString());
        }

    }
}
