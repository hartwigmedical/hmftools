package com.hartwig.hmftools.isofox.neo;

import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWN;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_PAIR;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UP;
import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.fusion.FusionCommon.switchStream;
import static com.hartwig.hmftools.common.codon.Nucleotides.reverseComplementBases;
import static com.hartwig.hmftools.common.bam.CigarUtils.cigarFromStr;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.generateMappedCoords;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.isofox.neo.NeoFragmentSupport.EXACT_MATCH;

import java.util.List;

import com.hartwig.hmftools.common.neo.NeoEpitopeFile;
import com.hartwig.hmftools.common.neo.RnaNeoEpitope;

import org.apache.commons.compress.utils.Lists;

import htsjdk.samtools.Cigar;

public class NeoEpitopeData
{
    public final NeoEpitopeFile Source;

    public final String[] Chromosomes;
    public final int[] Positions;
    public final byte[] Orientations;

    public final List<String>[] Transcripts;
    public final List<int[]>[] CodingBaseCoords;

    private final NeoFragmentSupport mFragmentSupport;

    public NeoEpitopeData(final NeoEpitopeFile source)
    {
        Source = source;

        Chromosomes = new String[FS_PAIR];
        Positions = new int[FS_PAIR];
        Orientations = new byte[FS_PAIR];
        Source.extractLocationData(Chromosomes, Positions, Orientations);

        Transcripts = new List[FS_PAIR];
        Transcripts[FS_UP] = Lists.newArrayList();
        Transcripts[FS_DOWN] = Lists.newArrayList();
        Source.extractTranscriptNames(Transcripts[FS_UP], Transcripts[FS_DOWN]);

        CodingBaseCoords = new List[FS_PAIR];

        for(int fs = FS_UP; fs <= FS_DOWN; ++fs)
        {
            CodingBaseCoords[fs] = Lists.newArrayList();

            if(Source.CodingBaseCigars[fs].isEmpty())
                continue;

            final Cigar cigar = cigarFromStr(Source.CodingBaseCigars[fs]);
            CodingBaseCoords[fs].addAll(generateMappedCoords(cigar, Source.CodingBasePositions[fs][SE_START]));
        }

        mFragmentSupport = new NeoFragmentSupport();
    }

    public boolean isFusion() { return Source.VariantType.isFusion(); }
    public boolean isPointMutation() { return Source.VariantType.isPointMutation(); }

    public boolean singleGene() { return Source.GeneIds[FS_UP].equals(Source.GeneIds[FS_DOWN]); }
    public boolean singleChromosome() { return Chromosomes[FS_UP].equals(Chromosomes[FS_DOWN]); }
    public boolean posStrand(int stream) { return (stream == FS_UP) == (Orientations[stream] == POS_ORIENT); }

    public NeoFragmentSupport getFragmentSupport() { return mFragmentSupport; }

    public void setOrientation(final byte geneStrand)
    {
        if(geneStrand == POS_STRAND)
        {
            Orientations[FS_UP] = POS_ORIENT;
            Orientations[FS_DOWN] = NEG_ORIENT;
        }
        else
        {
            Orientations[FS_UP] = NEG_ORIENT;
            Orientations[FS_DOWN] = POS_ORIENT;
        }
    }

    public int[] getCodingBaseRange(int stream) { return Source.CodingBasePositions[stream]; }

    public String getFullCodingBases(int streamPerspective)
    {
        if(!isFusion() || (Orientations[FS_UP] != Orientations[FS_DOWN] && singleChromosome()))
        {
            if(Orientations[FS_UP] == POS_ORIENT)
                return Source.CodingBases[FS_UP] + Source.CodingBases[FS_DOWN];
            else
                return Source.CodingBases[FS_DOWN] + Source.CodingBases[FS_UP];
        }

        if(streamPerspective == FS_UP)
        {
            if(Orientations[FS_UP] == POS_ORIENT)
                return Source.CodingBases[FS_UP] + reverseComplementBases(Source.CodingBases[FS_DOWN]);
            else
                return reverseComplementBases(Source.CodingBases[FS_DOWN]) + Source.CodingBases[FS_UP];
        }
        else
        {
            if(Orientations[FS_DOWN] == POS_ORIENT)
                return Source.CodingBases[FS_DOWN] + reverseComplementBases(Source.CodingBases[FS_UP]);
            else
                return reverseComplementBases(Source.CodingBases[FS_UP]) + Source.CodingBases[FS_DOWN];
        }
    }

    public boolean isDeletionFusion()
    {
        if(!isFusion())
            return false;

        if(!singleChromosome())
            return false;

        if(Positions[FS_UP] < Positions[FS_DOWN])
        {
            return Orientations[FS_UP] == POS_ORIENT && Orientations[FS_DOWN] == NEG_ORIENT;
        }
        else
        {
            return Orientations[FS_UP] == NEG_ORIENT && Orientations[FS_DOWN] == POS_ORIENT;
        }
    }

    public String getFusionSoftClippedBases(int streamPerspective)
    {
        if(Orientations[FS_UP] != Orientations[FS_DOWN])
        {
            return Source.CodingBases[switchStream(streamPerspective)];
        }
        else
        {
            return reverseComplementBases(Source.CodingBases[switchStream(streamPerspective)]);
        }
    }

    public RnaNeoEpitope asRnaFile()
    {
        return new RnaNeoEpitope(
                Source.Id, Source.VariantType, Source.VariantInfo, mFragmentSupport.NovelFragments[EXACT_MATCH],
                mFragmentSupport.RefBaseDepth[FS_UP], mFragmentSupport.RefBaseDepth[FS_DOWN]);
    }
}
