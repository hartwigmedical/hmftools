package com.hartwig.hmftools.linx.neoepitope;

import static com.hartwig.hmftools.common.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.linx.neoepitope.AminoAcidConverter.STOP_SYMBOL;
import static com.hartwig.hmftools.linx.neoepitope.AminoAcidConverter.convertDnaCodonToAminoAcid;
import static com.hartwig.hmftools.linx.neoepitope.AminoAcidConverter.reverseStrandBases;
import static com.hartwig.hmftools.linx.neoepitope.AminoAcidConverter.isStopCodon;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_END;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_START;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.structural.annotation.ExonData;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneAnnotation;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneFusion;
import com.hartwig.hmftools.common.variant.structural.annotation.Transcript;
import com.hartwig.hmftools.common.variant.structural.annotation.TranscriptData;
import com.hartwig.hmftools.linx.gene.SvGeneTranscriptCollection;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class NeoEpitopeFinder
{
    private final SvGeneTranscriptCollection mGeneTransCache;
    private final RefGenomeInterface mRefGenome;
    private final String mOutputDir;
    private BufferedWriter mFileWriter;

    private final List<NeoEpitopeData> mNeoEpitopeResults;

    private static final int AMINO_ACID_REF_COUNT = 10;

    private static final Logger LOGGER = LogManager.getLogger(NeoEpitopeFinder.class);

    public NeoEpitopeFinder(final RefGenomeInterface refGenome, final SvGeneTranscriptCollection geneTransCache, final String outputDir)
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
            final Transcript upTrans = fusion.upstreamTrans();
            final Transcript downTrans = fusion.downstreamTrans();

            if(!upTrans.isCanonical() || !downTrans.isCanonical())
                continue;

            if(upTrans.nonCoding() || downTrans.nonCoding())
                continue;

            if(upTrans.preCoding() || downTrans.preCoding() || downTrans.isPromoter())
                continue;

            if(upTrans.isExonic() && downTrans.isExonic())
                continue;

            boolean isPhased = fusion.phaseMatched()
                    && fusion.getExonsSkipped(true) == 0 && fusion.getExonsSkipped(false) == 0;

            int upstreamPhaseOffset = upTrans.ExonUpstreamPhase;
            int downstreamPhaseOffset = upstreamPhaseOffset == 0 || !isPhased ? 0 : 3 - upstreamPhaseOffset;

            LOGGER.debug("fusion({}) phased({}) upPhaseOffset({}) downPhaseOffset({})",
                    fusion.name(), isPhased, upstreamPhaseOffset, downstreamPhaseOffset);

            String upstreamBases = getBaseString(upTrans, false, upstreamPhaseOffset);
            String downstreamBases = getBaseString(downTrans, !isPhased, downstreamPhaseOffset);

            // upstream strand 1, bases will be retreived from left to right (lower to higher), no need for any conversion
            // downstream strand 1, bases will be retreived from left to right (lower to higher), no need for any conversion
            // upstream strand -1, bases will be retreived from left to right (lower to higher), need to reverse and convert
            // downstream strand -1, bases will be retreived from left to right (lower to higher), need to reverse and convert

            int upStrand = upTrans.parent().Strand;
            int downStrand = downTrans.parent().Strand;

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

            LOGGER.debug("fusion({}) upBases({}) novelCodon({}) downBases({})",
                    fusion.name(), upstreamBases, checkTrimBases(novelCodonBases), checkTrimBases(downstreamBases));

            if(upstreamBases.isEmpty() || (downstreamBases.isEmpty() && novelCodonBases.isEmpty()))
                continue;

            final String upstreamRefAminoAcids = getAminoAcids(upstreamBases, false);
            final String novelAminoAcids = getAminoAcids(novelCodonBases, !isPhased);
            String downstreamRefAminoAcids = getAminoAcids(downstreamBases, !isPhased);

            LOGGER.debug("fusion({}) upAA({}) novel({}) downAA({})",
                    fusion.name(), upstreamRefAminoAcids, checkTrimBases(novelAminoAcids), checkTrimBases(downstreamRefAminoAcids));

            if(novelAminoAcids.equals(STOP_SYMBOL))
                downstreamRefAminoAcids = "";

            NeoEpitopeData neoEpData = ImmutableNeoEpitopeData.builder()
                    .fusion(fusion)
                    .upstreamAcids(upstreamRefAminoAcids)
                    .downstreamAcids(downstreamRefAminoAcids)
                    .novelAcid(novelAminoAcids)
                    .build();

            mNeoEpitopeResults.add(neoEpData);

            writeData(sampleId, neoEpData, fusion);
        }
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

    private String getBaseString(final Transcript transcript, boolean collectAllBases, int phaseOffset)
    {
        if(transcript.nonCoding())
            return "";

        final GeneAnnotation gene = transcript.parent();

        long breakPosition = gene.position();

        final TranscriptData transData = mGeneTransCache.getTranscriptData(gene.StableId, transcript.StableId);

        if(transData == null)
        {
            LOGGER.error("gene({}) transcript({}) data not found", gene.GeneName, transcript.StableId);
            return "";
        }

        int requiredBases = AMINO_ACID_REF_COUNT * 3 + phaseOffset;

        long codingStart = transcript.CodingStart;
        long codingEnd = transcript.CodingEnd;

        final List<ExonData> exonDataList = transData.exons();

        String baseString = "";

        if(gene.orientation() == -1)
        {
            // walk forwards through the exons, collecting up the required positions and bases
            for (int i = 0; i < exonDataList.size(); ++i)
            {
                final ExonData exon = exonDataList.get(i);

                if (exon.ExonStart < breakPosition)
                    continue;

                long posStart, posEnd;

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

                // stop at end of coding region
                if(posStart < codingStart)
                {
                    posStart = codingStart;
                    requiredBases = 0;
                }

                if(posEnd > codingEnd)
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

                if(exon.ExonEnd > breakPosition)
                    continue;

                long posStart, posEnd;

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

                if(posStart < codingStart)
                {
                    posStart = codingStart;
                    requiredBases = 0;
                }

                if(posEnd > codingEnd)
                {
                    posEnd = codingEnd;
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

    private void writeData(final String sampleId, final NeoEpitopeData data, final GeneFusion fusion)
    {
        if(mOutputDir.isEmpty())
            return;

        try
        {
            if(mFileWriter == null)
            {
                String outputFileName = mOutputDir + "SVA_NEO_EPITOPES.csv";

                mFileWriter = createBufferedWriter(outputFileName, false);

                mFileWriter.write("SampleId,Fusion,PhaseMatched");
                mFileWriter.write(",UpstreamAminoAcids,DownstreamAminoAcids,NovelAminoAcid");

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
                    sampleId, fusion.name(), fusion.phaseMatched()));

            mFileWriter.write(String.format(",%s,%s,%s",
                    data.upstreamAcids(), data.downstreamAcids(), data.novelAcid()));

            for(int se = SE_START; se <= SE_END; ++se)
            {
                boolean isUpstream = (se == SE_START);
                final Transcript trans = isUpstream ? fusion.upstreamTrans() : fusion.downstreamTrans();
                final GeneAnnotation gene = trans.parent();

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
            LOGGER.error("error writing kataegis output file: {}", e.toString());
        }
    }

    public void close()
    {
        closeBufferedWriter(mFileWriter);
    }

}
