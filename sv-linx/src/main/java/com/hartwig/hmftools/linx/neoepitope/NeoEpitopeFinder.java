package com.hartwig.hmftools.linx.neoepitope;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.linx.LinxConfig.LNX_LOGGER;
import static com.hartwig.hmftools.linx.fusion.FusionFinder.isIrrelevantSameGene;
import static com.hartwig.hmftools.linx.fusion.FusionFinder.validFusionTranscript;
import static com.hartwig.hmftools.linx.neoepitope.AminoAcidConverter.STOP_SYMBOL;
import static com.hartwig.hmftools.linx.neoepitope.AminoAcidConverter.convertDnaCodonToAminoAcid;
import static com.hartwig.hmftools.linx.neoepitope.AminoAcidConverter.isStopCodon;
import static com.hartwig.hmftools.linx.neoepitope.AminoAcidConverter.reverseStrandBases;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.ensemblcache.ExonData;
import com.hartwig.hmftools.common.fusion.GeneAnnotation;
import com.hartwig.hmftools.common.fusion.GeneFusion;
import com.hartwig.hmftools.common.fusion.Transcript;
import com.hartwig.hmftools.common.ensemblcache.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;

public class NeoEpitopeFinder
{
    private final EnsemblDataCache mGeneTransCache;
    private final RefGenomeInterface mRefGenome;
    private final String mOutputDir;
    private BufferedWriter mFileWriter;

    private final List<NeoEpitopeData> mNeoEpitopeResults;

    private static final int AMINO_ACID_REF_COUNT = 10;

    public NeoEpitopeFinder(final RefGenomeInterface refGenome, final EnsemblDataCache geneTransCache, final String outputDir)
    {
        mGeneTransCache = geneTransCache;
        mRefGenome = refGenome;
        mNeoEpitopeResults = Lists.newArrayList();

        mOutputDir = outputDir;
        mFileWriter = null;
    }

    public final List<NeoEpitopeData> getResults() { return mNeoEpitopeResults; }

    public void reportNeoEpitopes(final String sampleId, final List<GeneFusion> fusions)
    {
        mNeoEpitopeResults.clear();

        for(final GeneFusion fusion : fusions)
        {
            if(isDuplicate(fusion))
                continue;

            final Transcript upTrans = fusion.upstreamTrans();
            final Transcript downTrans = fusion.downstreamTrans();

            if(!upTrans.isCanonical() || !downTrans.isCanonical())
                continue;

            if(!validTranscripts(upTrans, downTrans))
                continue;

            boolean isPhased = fusion.phaseMatched()
                    && fusion.getExonsSkipped(true) == 0 && fusion.getExonsSkipped(false) == 0;

            int upstreamPhaseOffset = upTrans.ExonUpstreamPhase;
            int downstreamPhaseOffset = upstreamPhaseOffset == 0 || !isPhased ? 0 : 3 - upstreamPhaseOffset;

            LNX_LOGGER.debug("fusion({}) SVs({} & {}) phased({}) upPhaseOffset({}) downPhaseOffset({})",
                    fusion.name(), fusion.svId(true), fusion.svId(false),
                    isPhased, upstreamPhaseOffset, downstreamPhaseOffset);


            String upstreamBases = getBaseString(upTrans, getTranscriptData(upTrans), false, upstreamPhaseOffset);
            final TranscriptData downTransData = getTranscriptData(downTrans);
            String downstreamBases = getBaseString(downTrans, downTransData, !isPhased, downstreamPhaseOffset);
            int nmdBaseCount = calcNonMediatedDecayBases(downTrans.gene(), downTransData);

            // upstream strand 1, bases will be retreived from left to right (lower to higher), no need for any conversion
            // downstream strand 1, bases will be retreived from left to right (lower to higher), no need for any conversion
            // upstream strand -1, bases will be retreived from left to right (lower to higher), need to reverse and convert
            // downstream strand -1, bases will be retreived from left to right (lower to higher), need to reverse and convert

            int upStrand = upTrans.gene().Strand;
            int downStrand = downTrans.gene().Strand;

            // correct for strand
            if(upStrand == -1)
                upstreamBases = reverseStrandBases(upstreamBases);

            if(downStrand == -1)
                downstreamBases = reverseStrandBases(downstreamBases);

            String novelCodonBases = "";

            if(upstreamPhaseOffset > upstreamBases.length() || downstreamPhaseOffset > downstreamBases.length())
                continue;

            // if upstream ends on a phase other than 0, need to take the bases from the downstream gene to make a novel codon

            if(upstreamPhaseOffset > 0)
            {
                // take the last 1 or 2 bases from the end of upstream gene's section
                novelCodonBases = upstreamBases.substring(upstreamBases.length() - upstreamPhaseOffset);
                upstreamBases = upstreamBases.substring(0, upstreamBases.length() - upstreamPhaseOffset);
            }

            if(isPhased)
            {
                novelCodonBases += downstreamBases.substring(0, downstreamPhaseOffset);
                downstreamBases = downstreamBases.substring(downstreamPhaseOffset);
            }
            else
            {
                novelCodonBases += downstreamBases;
                downstreamBases = "";
            }

            LNX_LOGGER.debug("fusion({}) upBases({}) novelCodon({}) downBases({}) downNmdBases({})",
                    fusion.name(), upstreamBases, checkTrimBases(novelCodonBases), checkTrimBases(downstreamBases), nmdBaseCount);

            if(upstreamBases.isEmpty() || (downstreamBases.isEmpty() && novelCodonBases.isEmpty()))
                continue;

            final String upstreamRefAminoAcids = getAminoAcids(upstreamBases, false);
            final String novelAminoAcids = getAminoAcids(novelCodonBases, !isPhased);
            String downstreamRefAminoAcids = getAminoAcids(downstreamBases, !isPhased);

            LNX_LOGGER.debug("fusion({}) upAA({}) novel({}) downAA({})",
                    fusion.name(), upstreamRefAminoAcids, checkTrimBases(novelAminoAcids), checkTrimBases(downstreamRefAminoAcids));

            if(novelAminoAcids.equals(STOP_SYMBOL))
                downstreamRefAminoAcids = "";

            NeoEpitopeData neoEpData = ImmutableNeoEpitopeData.builder()
                    .fusion(fusion)
                    .upstreamAcids(upstreamRefAminoAcids)
                    .downstreamAcids(downstreamRefAminoAcids)
                    .novelAcid(novelAminoAcids)
                    .downstreamNmdBases(nmdBaseCount)
                    .build();

            mNeoEpitopeResults.add(neoEpData);

            writeData(sampleId, neoEpData, fusion);
        }
    }

    private boolean isDuplicate(final GeneFusion fusion)
    {
        return mNeoEpitopeResults.stream()
                .anyMatch(x -> x.fusion().svId(true) == fusion.svId(true)
                && x.fusion().svId(false) == fusion.svId(false));
    }

    private static String checkTrimBases(final String bases)
    {
        if(bases.length() < 50)
            return bases;

        return bases.substring(0, 50) + "...";
    }

    private String getAminoAcids(final String baseString, boolean checkStopCodon)
    {
        if(baseString.length() < 3)
            return "";

        String aminoAcidStr = "";
        int index = 0;
        while(index <= baseString.length() - 3)
        {
            String codonBases = baseString.substring(index, index + 3);

            if(isStopCodon(codonBases) && checkStopCodon)
                break;

            String aminoAcid = convertDnaCodonToAminoAcid(codonBases);

            aminoAcidStr += aminoAcid;
            index += 3;
        }

        return aminoAcidStr;
    }

    private TranscriptData getTranscriptData(final Transcript transcript)
    {
        final TranscriptData transData = mGeneTransCache.getTranscriptData(transcript.gene().StableId, transcript.StableId);

        if(transData == null)
        {
            LNX_LOGGER.error("gene({}) transcript({}) data not found", transcript.gene().GeneName, transcript.StableId);
            return null;
        }

        return transData;
    }

    private String getBaseString(final Transcript transcript, final TranscriptData transData, boolean collectAllBases, int phaseOffset)
    {
        if(transcript.nonCoding())
            return "";

        final GeneAnnotation gene = transcript.gene();
        int breakPosition = gene.position();

        int requiredBases = AMINO_ACID_REF_COUNT * 3 + phaseOffset;

        int codingStart = transcript.CodingStart;
        int codingEnd = transcript.CodingEnd;
        boolean postCoding = transcript.postCoding();

        final List<ExonData> exonDataList = transData.exons();

        String baseString = "";

        if(gene.orientation() == -1)
        {
            // walk forwards through the exons, collecting up the required positions and bases
            for (int i = 0; i < exonDataList.size(); ++i)
            {
                final ExonData exon = exonDataList.get(i);

                if (exon.ExonStart < breakPosition || (!postCoding && exon.ExonStart < codingStart))
                    continue;

                int posStart, posEnd;

                if (collectAllBases || requiredBases > exon.ExonEnd - exon.ExonStart + 1)
                {
                    posStart = exon.ExonStart;
                    posEnd = exon.ExonEnd;
                    requiredBases -= (exon.ExonEnd - exon.ExonStart + 1);
                }
                else
                {
                    posStart = exon.ExonStart;
                    posEnd = exon.ExonStart + requiredBases - 1;
                    requiredBases = 0;
                }

                // stop at end of coding region unless the breakend started past it
                if(!postCoding && posEnd > codingEnd)
                {
                    posEnd = codingEnd;
                    requiredBases = 0;
                }

                if(posEnd < posStart)
                    continue;

                baseString += mRefGenome.getBaseString(gene.chromosome(), posStart, posEnd);

                if (requiredBases <= 0)
                    break;
            }
        }
        else
        {
            for(int i = exonDataList.size() - 1; i >= 0; --i)
            {
                final ExonData exon = exonDataList.get(i);

                if(exon.ExonEnd > breakPosition || (!postCoding && exon.ExonEnd > codingEnd))
                    continue;

                int posStart, posEnd;

                if(collectAllBases || requiredBases > exon.ExonEnd - exon.ExonStart + 1)
                {
                    posEnd = exon.ExonEnd;
                    posStart = exon.ExonStart;
                    requiredBases -= (exon.ExonEnd - exon.ExonStart + 1);
                }
                else
                {
                    posEnd = exon.ExonEnd;
                    posStart = exon.ExonEnd - requiredBases + 1;
                    requiredBases = 0;
                }

                if(!postCoding && posStart < codingStart)
                {
                    posStart = codingStart;
                    requiredBases = 0;
                }

                if(posEnd < posStart)
                    continue;

                // add in reverse since walking backwards through the exons
                baseString = mRefGenome.getBaseString(gene.chromosome(), posStart, posEnd) + baseString;

                if(requiredBases <= 0)
                    break;
            }

        }

        return baseString;
    }

    private int calcNonMediatedDecayBases(final GeneAnnotation gene, final TranscriptData transData)
    {
        final List<ExonData> exonDataList = transData.exons();

        int breakPosition = gene.position();

        int exonicBaseCount = 0;

        if(gene.orientation() == -1)
        {
            // walk forwards through the exons, collecting up the required positions and bases
            for (int i = 0; i < exonDataList.size(); ++i)
            {
                final ExonData exon = exonDataList.get(i);

                if(i == exonDataList.size() - 1)
                    break;

                if (breakPosition > exon.ExonEnd)
                    continue;

                exonicBaseCount += exon.ExonEnd - max(breakPosition, exon.ExonStart) + 1;
            }
        }
        else
        {
            for(int i = exonDataList.size() - 1; i >= 0; --i)
            {
                final ExonData exon = exonDataList.get(i);

                if(i == 0)
                    break;

                if(breakPosition < exon.ExonStart)
                    continue;

                exonicBaseCount += min(breakPosition, exon.ExonEnd) - exon.ExonStart + 1;
            }
        }

        return exonicBaseCount;
    }

    public void checkFusions(List<GeneFusion> existingFusions, List<GeneAnnotation> genesList1,List<GeneAnnotation> genesList2)
    {
        for (final GeneAnnotation gene1 : genesList1)
        {
            boolean startUpstream = gene1.isUpstream();

            final Transcript trans1 = gene1.transcripts().stream().filter(Transcript::isCanonical).findFirst().orElse(null);

            if (trans1 == null)
                continue;

            for (final GeneAnnotation gene2 : genesList2)
            {
                boolean endUpstream = gene2.isUpstream();

                if (startUpstream == endUpstream)
                    continue;

                final Transcript trans2 = gene2.transcripts().stream().filter(Transcript::isCanonical).findFirst().orElse(null);

                if (trans2 == null)
                    continue;

                final Transcript upstreamTrans = startUpstream ? trans1 : trans2;
                final Transcript downstreamTrans = upstreamTrans == trans1 ? trans2 : trans1;

                if(existingFusions.stream().anyMatch(x -> x.upstreamTrans() == upstreamTrans && x.downstreamTrans() == downstreamTrans))
                    continue;

                GeneFusion newFusion = checkNeoEpitopeFusion(upstreamTrans, downstreamTrans);

                if(newFusion != null)
                    existingFusions.add(newFusion);
            }
        }
    }

    private boolean validTranscripts(final Transcript upstreamTrans, final Transcript downstreamTrans)
    {
        if(!validFusionTranscript(upstreamTrans))
            return false;

        if(downstreamTrans.nonCoding() || downstreamTrans.ExonMax == 1)
            return false;

        if(upstreamTrans.preCoding() || upstreamTrans.ExonUpstreamPhase == -1)
            return false;

        if(upstreamTrans.isExonic() && downstreamTrans.isExonic()) // rare and too complicated for now
            return false;

        return true;
    }

    private GeneFusion checkNeoEpitopeFusion(final Transcript upstreamTrans, final Transcript downstreamTrans)
    {
        if(!validTranscripts(upstreamTrans, downstreamTrans))
            return null;

        if (isIrrelevantSameGene(upstreamTrans, downstreamTrans))
            return null;

        boolean phaseMatched = upstreamTrans.ExonUpstreamPhase == downstreamTrans.ExonDownstreamPhase;

        GeneFusion fusion = new GeneFusion(upstreamTrans, downstreamTrans, phaseMatched);
        fusion.setNeoEpitopeOnly(true);
        return fusion;
    }

    private void writeData(final String sampleId, final NeoEpitopeData data, final GeneFusion fusion)
    {
        if(mOutputDir.isEmpty())
            return;

        try
        {
            if(mFileWriter == null)
            {
                String outputFileName = mOutputDir + "LNX_NEO_EPITOPES.csv";

                mFileWriter = createBufferedWriter(outputFileName, false);

                mFileWriter.write("SampleId,Fusion,SameGene");
                mFileWriter.write(",UpstreamAminoAcids,DownstreamAminoAcids,NovelAminoAcid,NMDBases");

                for(int se = SE_START; se <= SE_END; ++se)
                {
                    String upDown = se == SE_START ? "Up" : "Down";

                    String fieldsStr = ",SvId" + upDown;
                    fieldsStr += ",Chr" + upDown;
                    fieldsStr += ",Pos" + upDown;
                    fieldsStr += ",Orient" + upDown;
                    fieldsStr += ",Trans" + upDown;
                    fieldsStr += ",Strand" + upDown;
                    fieldsStr += ",RegionType" + upDown;
                    fieldsStr += ",CodingType" + upDown;
                    fieldsStr += ",Exon" + upDown;
                    fieldsStr += ",Phase" + upDown;
                    mFileWriter.write(fieldsStr);
                }

                mFileWriter.newLine();
            }

            mFileWriter.write(String.format("%s,%s,%s",
                    sampleId, fusion.name(), fusion.upstreamTrans().geneName().equals(fusion.downstreamTrans().geneName())));

            mFileWriter.write(String.format(",%s,%s,%s,%d",
                    data.upstreamAcids(), data.downstreamAcids(), data.novelAcid(), data.downstreamNmdBases()));

            for(int se = SE_START; se <= SE_END; ++se)
            {
                boolean isUpstream = (se == SE_START);
                final Transcript trans = isUpstream ? fusion.upstreamTrans() : fusion.downstreamTrans();
                final GeneAnnotation gene = trans.gene();

                mFileWriter.write(String.format(",%d,%s,%d,%d",
                        gene.id(), gene.chromosome(), gene.position(), gene.orientation()));

                mFileWriter.write(String.format(",%s,%d,%s,%s,%d,%d",
                        trans.StableId, gene.Strand, trans.regionType(), trans.codingType(),
                        trans.nextSpliceExonRank(), trans.nextSpliceExonPhase()));
            }

            mFileWriter.newLine();
        }
        catch (final IOException e)
        {
            LNX_LOGGER.error("error writing kataegis output file: {}", e.toString());
        }
    }

    public void close()
    {
        closeBufferedWriter(mFileWriter);
    }

}
