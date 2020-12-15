package com.hartwig.hmftools.svtools.neo;

import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWNSTREAM;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UPSTREAM;
import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.fusion.TranscriptRegionType.EXONIC;
import static com.hartwig.hmftools.svtools.neo.NeoConfig.AMINO_ACID_REF_COUNT;
import static com.hartwig.hmftools.svtools.neo.NeoEpitopeAnnotator.IM_LOGGER;
import static com.hartwig.hmftools.svtools.neo.NeoUtils.getCodingBases;
import static com.hartwig.hmftools.svtools.neo.NeoUtils.setTranscriptCodingData;
import static com.hartwig.hmftools.svtools.neo.NeoUtils.setTranscriptContext;

import com.hartwig.hmftools.common.ensemblcache.TranscriptData;
import com.hartwig.hmftools.common.fusion.TranscriptRegionType;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.neo.NeoEpitopeFusion;

public class SvNeoEpitope extends NeoEpitope
{
    private final NeoEpitopeFusion mSvFusion;

    public SvNeoEpitope(final NeoEpitopeFusion fusion)
    {
        super();
        mSvFusion = fusion;
    }

    public int position(int stream)
    {
        return mSvFusion.Positions[stream];
    }

    public String chromosome(int stream)
    {
        return mSvFusion.Chromosomes[stream];
    }

    public void setTranscriptData(final TranscriptData upTransData, final TranscriptData downTransData)
    {
        TransData[FS_UPSTREAM] = upTransData;
        TransData[FS_DOWNSTREAM] = downTransData;

        for(int fs = FS_UPSTREAM; fs <= FS_DOWNSTREAM; ++fs)
        {
            setTranscriptContext(this, TransData[fs], position(fs), fs);
            int insSeqLength = fs == FS_UPSTREAM ? mSvFusion.InsertSequence.length() : 0;
            setTranscriptCodingData(this, TransData[fs], position(fs), insSeqLength, fs);
        }
    }

    public void extractCodingBases(final RefGenomeInterface refGenome, int requiredAminoAcids)
    {
        // get the number of bases from the fusion junction as required by the required amino acid count

        // for fusions, the upstream bases include any inserted bases if the region is exonic
        // phasing already takes insert sequence into account, so just adjust for the number of bases required

        int upPhaseOffset = Phases[FS_UPSTREAM];
        int downPhaseOffset = getDownstreamPhaseOffset();

        int upRequiredBases = requiredAminoAcids * 3 + upPhaseOffset;
        String codingInsSequence = "";

        if(RegionType[FS_UPSTREAM] == TranscriptRegionType.EXONIC && !mSvFusion.InsertSequence.isEmpty())
        {
            upRequiredBases -= mSvFusion.InsertSequence.length();
            codingInsSequence = mSvFusion.InsertSequence;
        }

        CodingBases[FS_UPSTREAM] = getCodingBases(
                refGenome, TransData[FS_UPSTREAM], chromosome(FS_UPSTREAM), position(FS_UPSTREAM), orientation(FS_UPSTREAM),
                upRequiredBases, true);

        if(!codingInsSequence.isEmpty())
        {
            if(TransData[FS_UPSTREAM].Strand == POS_STRAND)
                CodingBases[FS_UPSTREAM] += codingInsSequence;
            else
                CodingBases[FS_UPSTREAM]= codingInsSequence + CodingBases[FS_UPSTREAM];
        }

        boolean canStartInExon = RegionType[FS_UPSTREAM] == TranscriptRegionType.EXONIC;
        int downRequiredBases = requiredAminoAcids * 3 + downPhaseOffset;

        CodingBases[FS_DOWNSTREAM] = getCodingBases(
                refGenome, TransData[FS_DOWNSTREAM], chromosome(FS_DOWNSTREAM), position(FS_DOWNSTREAM), orientation(FS_DOWNSTREAM),
                downRequiredBases, canStartInExon);

        IM_LOGGER.debug("ne({}) phased({} up={} down={}) reqBases(up={} down={}) insSeqLen({})",
                toString(), phaseMatched(), upPhaseOffset, downPhaseOffset, upRequiredBases, downRequiredBases, codingInsSequence.length());
    }

    public String toString()
    {
        return String.format("fusion up(%s: %s:%d:%d) down(%s: %s:%d:%d)",
                mSvFusion.GeneNames[FS_UPSTREAM], mSvFusion.Chromosomes[FS_UPSTREAM],
                mSvFusion.Positions[FS_UPSTREAM], mSvFusion.Orientations[FS_UPSTREAM],
                mSvFusion.GeneNames[FS_DOWNSTREAM], mSvFusion.Chromosomes[FS_DOWNSTREAM],
                mSvFusion.Positions[FS_DOWNSTREAM], mSvFusion.Orientations[FS_DOWNSTREAM]);
    }


}
