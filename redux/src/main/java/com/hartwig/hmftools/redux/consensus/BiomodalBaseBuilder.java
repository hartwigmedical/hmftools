package com.hartwig.hmftools.redux.consensus;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.BASE_MODIFICATIONS_ATTRIBUTE;
import static com.hartwig.hmftools.common.codon.Nucleotides.swapDnaBase;
import static com.hartwig.hmftools.common.sequencing.BiomodalBamUtils.MODC_ANNOTATION;
import static com.hartwig.hmftools.common.sequencing.BiomodalBamUtils.getMMValueFromModCReadIndices;
import static com.hartwig.hmftools.common.sequencing.BiomodalBamUtils.getModCReadIndices;
import static com.hartwig.hmftools.redux.consensus.BaseBuilder.INVALID_POSITION;
import static com.hartwig.hmftools.redux.consensus.BaseBuilder.NO_BASE;

import static htsjdk.samtools.CigarOperator.I;
import static htsjdk.samtools.CigarOperator.M;

import java.util.Collection;
import java.util.List;
import java.util.SortedSet;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.redux.BaseQualAdjustment;

import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

public class BiomodalBaseBuilder extends NonStandardBaseBuilder
{
    private final BaseBuilder mBaseBuilder;

    public BiomodalBaseBuilder(final RefGenome refGenome)
    {
        super(refGenome);

        mBaseBuilder = new BaseBuilder(refGenome, null);
    }

    @Override
    public void buildConsensusRead(final List<SAMRecord> reads, final ConsensusState consensusState, final boolean hasIndels)
    {
        String chromosome = reads.get(0).getReferenceName();
        mChromosomeLength = mRefGenome.getChromosomeLength(chromosome);

        List<List<AnnotatedBase>> annotatedReads = annotateReads(reads);
        for(int i = 0; i < reads.size(); i++)
        {
            addModCAnnotation(reads.get(i), annotatedReads.get(i));
        }

        Collection<List<AnnotatedBase>> alignment = alignAnnotatedReads(annotatedReads);

        List<AnnotatedBase> consensusBases = getConsensusBases(
                alignment, records -> determineConsensus(chromosome, consensusState.IsForward, records));

        finalizeConsensusState(mRefGenome, reads, consensusState, hasIndels, consensusBases);

        SortedSet<Integer> modCReadIndices = Sets.newTreeSet();
        for(int i = 0; i < consensusBases.size(); i++)
        {
            if(consensusBases.get(i).Annotations.contains(MODC_ANNOTATION))
            {
                modCReadIndices.add(i);
            }
        }

        String mmValue = getMMValueFromModCReadIndices(consensusState.Bases, modCReadIndices, consensusState.IsForward);
        consensusState.Attributes.put(BASE_MODIFICATIONS_ATTRIBUTE, mmValue);
    }

    private AnnotatedBase determineConsensus(final String chromosome, boolean isForward, final List<AnnotatedBase> bases)
    {
        CigarOperator consensusOp = getConsensusCigarOp(bases);
        ExtendedRefPos pos = bases.get(0).Pos;
        byte[] locationBases = new byte[bases.size()];
        byte[] locationQuals = new byte[bases.size()];
        for(int i = 0; i < bases.size(); i++)
        {
            locationBases[i] = bases.get(i).Base;
            locationQuals[i] = bases.get(i).Qual;
        }

        int basePosition = pos.RefPos;
        if(consensusOp != M || basePosition < 1 || basePosition > mChromosomeLength)
        {
            basePosition = INVALID_POSITION;
        }

        byte[] consensusBaseAndQual = mBaseBuilder.determineBaseAndQual(locationBases, locationQuals, chromosome, basePosition);
        byte consensusBase = consensusBaseAndQual[0] == NO_BASE ? ANY_BASE : consensusBaseAndQual[0];
        byte consensusQual = BaseQualAdjustment.adjustBaseQual(consensusBaseAndQual[1]);

        if(consensusQual == (byte) 0 && consensusOp == I)
        {
            return null;
        }

        AnnotatedBase consensus = new AnnotatedBase(pos, consensusBase, consensusQual, consensusOp);

        byte targetBase = isForward ? (byte) 'C' : (byte) swapDnaBase('C');
        if(consensusBase != targetBase)
        {
            return consensus;
        }

        int nonModCQualTotal = 0;
        int modCQualTotal = 0;
        for(AnnotatedBase annotatedBase : bases)
        {
            if(annotatedBase.Base != targetBase)
            {
                continue;
            }

            if(annotatedBase.Annotations.contains(MODC_ANNOTATION))
            {
                modCQualTotal += annotatedBase.Qual;
                continue;
            }

            nonModCQualTotal += annotatedBase.Qual;
        }

        if(modCQualTotal <= nonModCQualTotal)
        {
            return consensus;
        }

        consensus.Annotations.add(MODC_ANNOTATION);
        return consensus;
    }

    private static void addModCAnnotation(final SAMRecord read, List<AnnotatedBase> annotatedRead)
    {
        SortedSet<Integer> modCReadIndices = getModCReadIndices(read);
        for(int readIndex : modCReadIndices)
        {
            annotatedRead.get(readIndex).Annotations.add(MODC_ANNOTATION);
        }
    }
}
