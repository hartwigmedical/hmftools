package com.hartwig.hmftools.neo.epitope;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.gene.CodingBaseData.PHASE_1;
import static com.hartwig.hmftools.common.gene.CodingBaseData.PHASE_2;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWN;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UP;
import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.gene.TranscriptRegionType.EXONIC;
import static com.hartwig.hmftools.common.gene.TranscriptRegionType.INTRONIC;
import static com.hartwig.hmftools.common.gene.TranscriptRegionType.UPSTREAM;
import static com.hartwig.hmftools.common.gene.TranscriptUtils.tickPhaseForward;
import static com.hartwig.hmftools.common.neo.NeoEpitopeFile.fusionInfo;
import static com.hartwig.hmftools.common.neo.NeoEpitopeType.INFRAME_FUSION;
import static com.hartwig.hmftools.common.neo.NeoEpitopeType.OUT_OF_FRAME_FUSION;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsWithin;
import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.genome.region.Orientation.ORIENT_REV;
import static com.hartwig.hmftools.common.genome.region.Orientation.ORIENT_FWD;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.gene.TranscriptRegionType;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.neo.NeoEpitopeFusion;
import com.hartwig.hmftools.common.neo.NeoEpitopeType;

import com.google.common.collect.Lists;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;

public class SvNeoEpitope extends NeoEpitope
{
    private final NeoEpitopeFusion mSvFusion;
    private final int[] mSkippedSpliceAcceptorDonors;

    public SvNeoEpitope(final NeoEpitopeFusion fusion)
    {
        super();
        mSvFusion = fusion;
        mSkippedSpliceAcceptorDonors = new int[] { 0, 0 };
    }

    public int position(int stream)
    {
        return mSvFusion.Positions[stream];
    }

    public String chromosome(int stream)
    {
        return mSvFusion.Chromosomes[stream];
    }
    public String geneName(int stream) { return mSvFusion.GeneNames[stream]; }

    public NeoEpitopeType variantType()
    {
        return phaseMatched() ? INFRAME_FUSION : OUT_OF_FRAME_FUSION;
    }

    public String variantInfo()
    {
        return fusionInfo(mSvFusion.Chromosomes, mSvFusion.Positions, mSvFusion.Orientations);
    }

    public double variantCopyNumber() { return mSvFusion.JunctionCopyNumber; }
    public double copyNumber() { return mSvFusion.CopyNumber; }
    public double subclonalLikelihood() { return 0; }

    public boolean phaseMatched()
    {
        if(RegionType[FS_UP] == EXONIC && RegionType[FS_DOWN] == EXONIC)
            return tickPhaseForward(Phases[FS_UP]) == Phases[FS_DOWN];
        else
            return Phases[FS_UP] == Phases[FS_DOWN];
    }

    public int unsplicedDistance()
    {
        int upEndCodingBase = mSvFusion.Orientations[FS_UP] == ORIENT_FWD ? ExtPositions[FS_UP][SE_END] : ExtPositions[FS_UP][SE_START];
        int upUnspliced = abs(position(FS_UP) - upEndCodingBase);
        int downEndCodingBase = mSvFusion.Orientations[FS_DOWN] == ORIENT_FWD ? ExtPositions[FS_DOWN][SE_END] : ExtPositions[FS_DOWN][SE_START];
        int downUnspliced = abs(position(FS_DOWN) - downEndCodingBase);
        return upUnspliced + downUnspliced + mSvFusion.ChainLength;
    }

    public int skippedDonors() { return mSkippedSpliceAcceptorDonors[FS_UP]; }
    public int skippedAcceptors() { return mSkippedSpliceAcceptorDonors[FS_DOWN]; }
    public String wildtypeAcids() { return ""; }

    public void setTranscriptData(final TranscriptData upTransData, final TranscriptData downTransData)
    {
        TransData[FS_UP] = upTransData;
        EpitopeUtils.setTranscriptContext(this, TransData[FS_UP], position(FS_UP), FS_UP);
        int insSeqLength = mSvFusion.InsertSequence.length();
        EpitopeUtils.setTranscriptCodingData(this, TransData[FS_UP], position(FS_UP), insSeqLength, FS_UP);

        TransData[FS_DOWN] = downTransData;

        // if the upstream context is intronic, then skip past any downstream exonic section

        EpitopeUtils.setTranscriptContext(this, TransData[FS_DOWN], position(FS_DOWN), FS_DOWN);
        EpitopeUtils.setTranscriptCodingData(this, TransData[FS_DOWN], position(FS_DOWN), 0, FS_DOWN);

        if(RegionType[FS_UP] == INTRONIC && RegionType[FS_DOWN] == EXONIC)
        {
            final ExonData exon = downTransData.exons().stream()
                    .filter(x -> positionWithin(position(FS_DOWN), x.Start, x.End))
                    .findFirst().orElse(null);

            if(exon != null)
            {
                Phases[FS_DOWN] = exon.PhaseEnd;
                ExonRank[FS_DOWN] = exon.Rank + 1;
            }
        }
        else if(RegionType[FS_DOWN] == UPSTREAM)
        {
            setDownstreamPhaseFromUtr();
        }
    }

    private int getUpstreamOpenCodonBases()
    {
        int phase = Phases[FS_UP];

        if(RegionType[FS_UP] == TranscriptRegionType.EXONIC && !mSvFusion.InsertSequence.isEmpty())
        {
            int insSeqLen = mSvFusion.InsertSequence.length();
            int preInsertPhase = (phase - insSeqLen) % 3;

            if(preInsertPhase >= 0)
                phase = preInsertPhase;
            else if(abs(preInsertPhase) == PHASE_1)
                phase = PHASE_2;
            else
                phase = PHASE_1;
        }

        return getUpstreamOpenCodonBases(phase);
    }

    private int getDownstreamOpenCodonBases()
    {
        int upOpenBases = Phases[FS_UP];
        return upOpenBases == 0 ? 0 : 3 - upOpenBases;
    }

    public void extractCodingBases(final RefGenomeInterface refGenome, int requiredAminoAcids)
    {
        // get the number of bases from the fusion junction as required by the required amino acid count

        // for fusions, the upstream bases include any inserted bases if the region is exonic
        // phasing already takes insert sequence into account, so just adjust for the number of bases required

        int upExtraBases = getUpstreamOpenCodonBases();
        int downExtraBases = getDownstreamOpenCodonBases();

        int upRequiredBases = requiredAminoAcids * 3 + upExtraBases;

        CodingBaseExcerpt cbExcerpt = EpitopeUtils.getUpstreamCodingBaseExcerpt(
                refGenome, TransData[FS_UP], chromosome(FS_UP), position(FS_UP), orientation(FS_UP), upRequiredBases);

        RawCodingBases[FS_UP] = cbExcerpt.Bases;

        ExtPositions[FS_UP][SE_START] = cbExcerpt.Positions[SE_START];
        ExtPositions[FS_UP][SE_END] = cbExcerpt.Positions[SE_END];

        List<CigarElement> upCigarElements = Lists.newArrayList();
        upCigarElements.addAll(cbExcerpt.CigarRef.getCigarElements());

        int codingInsSeqLen = 0;

        if(RegionType[FS_UP] == TranscriptRegionType.EXONIC && !mSvFusion.InsertSequence.isEmpty())
        {
            codingInsSeqLen = mSvFusion.InsertSequence.length();

            if(TransData[FS_UP].Strand == POS_STRAND)
            {
                upCigarElements.add(new CigarElement(codingInsSeqLen, CigarOperator.I));
                RawCodingBases[FS_UP] += mSvFusion.InsertSequence;
            }
            else
            {
                upCigarElements.add(0, new CigarElement(codingInsSeqLen, CigarOperator.I));
                RawCodingBases[FS_UP] = mSvFusion.InsertSequence + RawCodingBases[FS_UP];
            }
        }

        ExtCigars[FS_UP] = new Cigar();
        upCigarElements.forEach(x -> ExtCigars[FS_UP].add(x));
        ExtCodingBases[FS_UP] = RawCodingBases[FS_UP];

        NovelBaseIndex[FS_UP] = upExtraBases + codingInsSeqLen;
        NovelBaseIndex[FS_DOWN] = downExtraBases;

        boolean canStartInExon = RegionType[FS_UP] == TranscriptRegionType.EXONIC || upExtraBases > 0;
        int downRequiredBases = requiredAminoAcids * 3 + downExtraBases;

        cbExcerpt = EpitopeUtils.getDownstreamCodingBaseExcerpt(
                refGenome, TransData[FS_DOWN], chromosome(FS_DOWN), position(FS_DOWN), orientation(FS_DOWN),
                downRequiredBases, canStartInExon, true, !phaseMatched());

        if(cbExcerpt == null)
        {
            Valid = false;
            return;
        }

        RawCodingBases[FS_DOWN] = cbExcerpt.Bases;

        // call again to get restricted downstream bases
        if(!phaseMatched())
        {
            cbExcerpt = EpitopeUtils.getDownstreamCodingBaseExcerpt(
                    refGenome, TransData[FS_DOWN], chromosome(FS_DOWN), position(FS_DOWN), orientation(FS_DOWN),
                    downRequiredBases, canStartInExon, true, false);
        }

        if(cbExcerpt.Bases.isEmpty() || cbExcerpt.Positions[SE_START] == 0 || cbExcerpt.Positions[SE_END] == 0)
        {
            Valid = false;
            return;
        }

        ExtPositions[FS_DOWN][SE_START] = cbExcerpt.Positions[SE_START];
        ExtPositions[FS_DOWN][SE_END] = cbExcerpt.Positions[SE_END];
        ExtCodingBases[FS_DOWN] = cbExcerpt.Bases;
        ExtCigars[FS_DOWN] = cbExcerpt.CigarRef;

        NE_LOGGER.trace("ne({}) phased({} up={} down={}) reqBases(up={} down={}) insSeqLen({})",
                toString(), phaseMatched(), upExtraBases, downExtraBases, upRequiredBases, downRequiredBases, codingInsSeqLen);
    }

    public void checkStopLost(final RefGenomeInterface refGenome, int reqWildtypeAminoAcids) {}

    public void setSkippedSpliceSites(final EnsemblDataCache geneTransCache)
    {
        // between the fusion position and the nearest coding
        for(int fs = FS_UP; fs <= FS_DOWN; ++fs)
        {
            boolean findExonStart = mSvFusion.Orientations[fs] == ORIENT_REV;

            final int[] positionBoundaries = { 0, 0 };

            if(mSvFusion.Orientations[fs] == ORIENT_FWD)
            {
                positionBoundaries[SE_START] = ExtPositions[fs][SE_END];
                positionBoundaries[SE_END] = position(fs);
            }
            else
            {
                positionBoundaries[SE_START] = position(fs);
                positionBoundaries[SE_END] = ExtPositions[fs][SE_START];
            }

            final List<TranscriptData> transDataList = geneTransCache.getTranscripts(mSvFusion.GeneIds[fs]);
            final int stream = fs;

            final List<TranscriptData> candidateTransDataList = transDataList.stream()
                    .filter(x -> x.TransId != TransData[stream].TransId).collect(Collectors.toList());

            mSkippedSpliceAcceptorDonors[fs] = EpitopeUtils.findSkippedExonBoundaries(
                    candidateTransDataList, positionBoundaries, findExonStart, fs == FS_DOWN);
        }
    }

    private void setDownstreamPhaseFromUtr()
    {
        int breakpointPosition = position(FS_DOWN);
        TranscriptData transData = TransData[FS_DOWN];

        if(transData.posStrand())
        {
            for(ExonData exon : transData.exons())
            {
                if(exon.Rank == 1)
                    continue;

                if(exon.Start <= breakpointPosition)
                    continue;

                Phases[FS_DOWN] = exon.PhaseStart;
                ExonRank[FS_DOWN] = exon.Rank;
                break;
            }
        }
        else
        {
            for(int i = transData.exons().size() - 1; i >= 0; --i)
            {
                ExonData exon = transData.exons().get(i);

                if(exon.End >= breakpointPosition)
                    continue;

                Phases[FS_DOWN] = exon.PhaseEnd;
                ExonRank[FS_DOWN] = exon.Rank;
                break;
            }

        }
    }

    public String toString()
    {
        return String.format("fusion up(%s:%s %s:%d:%d) down(%s:%s %s:%d:%d) phased(%s) trans up(%s:%s) down(%s:%s)",
                mSvFusion.GeneNames[FS_UP], TransData[FS_UP].TransName, mSvFusion.Chromosomes[FS_UP],
                mSvFusion.Positions[FS_UP], mSvFusion.Orientations[FS_UP],
                mSvFusion.GeneNames[FS_DOWN], TransData[FS_DOWN].TransName, mSvFusion.Chromosomes[FS_DOWN],
                mSvFusion.Positions[FS_DOWN], mSvFusion.Orientations[FS_DOWN], phaseMatched(),
                CodingType[FS_UP], RegionType[FS_UP], CodingType[FS_DOWN], RegionType[FS_DOWN]);
    }

    public static boolean svIsNonDisruptiveInCodingTranscript(final int[] positions, final TranscriptData transcriptData)
    {
        if(transcriptData.nonCoding())
            return false;

        for(int i = 0; i < transcriptData.exons().size() - 1; ++i)
        {
            ExonData exon = transcriptData.exons().get(i);
            ExonData nextExon = transcriptData.exons().get(i + 1);

            if(positionsWithin(positions[SE_START], positions[SE_END], exon.End, nextExon.Start)
            && positionsWithin(exon.End, nextExon.Start, transcriptData.CodingStart, transcriptData.CodingEnd))
            {
                return true;
            }
        }

        return false;
    }
}
