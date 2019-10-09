package com.hartwig.hmftools.sage.context;

import java.util.function.Consumer;

import com.hartwig.hmftools.common.region.GenomeRegion;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class RefContextConsumer implements Consumer<SAMRecord> {

    private static final Logger LOGGER = LogManager.getLogger(RefContextConsumer.class);


    /*

    bcftools filter -i 'FORMAT/QUAL[1:0]>99 && FORMAT/AD[0:1] < 4' GIABvsSELFv004.sage.vcf.gz -O z -o GIABvsSELFv004.sage.filtered.vcf.gz
    bcftools filter -i 'FORMAT/QUAL[1:0]>99 && FORMAT/AD[0:1] < 4' COLO829v003.sage.vcf.gz -O z -o COLO829v003.sage.filtered.vcf.gz



    bgzip colo829.sage.vcf
    bcftools index colo829.sage.vcf.gz

    input=$1

    bcftools annotate -a /Users/jon/hmf/resources/GERMLINE_PON.vcf.gz -c GERMLINE_PON_COUNT ${input}.vcf.gz -O z -o ${input}.germline.vcf.gz
    bcftools index ${input}.germline.vcf.gz

    bcftools annotate -a /Users/jon/hmf/resources/SOMATIC_PON.vcf.gz -c SOMATIC_PON_COUNT ${input}.germline.vcf.gz -O z -o ${input}.pon.vcf.gz
    bcftools index ${input}.pon.vcf.gz

    bcftools annotate -a /Users/jon/hmf/resources/out_150_hg19.mappability.bed.gz -h /Users/jon/hmf/resources/mappability.hdr -c CHROM,FROM,TO,-,MAPPABILITY ${input}.pon.vcf.gz -O z -o ${input}.map.vcf.gz
    bcftools index ${input}.map.vcf.gz

    bcftools annotate -a all.somatic.snvs.vcf.gz -m PRE_STRELKA -c FILTER colo829.sage.map.vcf.gz -O z -o colo829.sage.pre.vcf.gz
    bcftools index colo829.sage.pre.vcf.gz
    bcftools annotate -a COLO829v003T.somatic_caller_post_processed.vcf.gz -m POST_STRELKA -c FILTER colo829.sage.pre.vcf.gz -O z -o colo829.sage.final.vcf.gz


    bcftools annotate -a all.somatic.snvs.vcf.gz -m PRE_STRELKA -c FILTER COLO829v003.sage.filtered.map.vcf.gz -O z -o COLO829v003.sage.filtered.pre.vcf.gz
    bcftools index COLO829v003.sage.filtered.pre.vcf.gz
    bcftools annotate -a COLO829v003T.somatic_caller_post_processed.vcf.gz -m POST_STRELKA -c FILTER COLO829v003.sage.filtered.pre.vcf.gz -O z -o COLO829v003.sage.filtered.final.vcf.gz


    */

    private static final String DISTANCE_FROM_REF_TAG = "NM";
    private final int minQuality;

    private final GenomeRegion bounds;
    private final IndexedFastaSequenceFile refGenome;
    private final RefContextCandidates candidates;
    private final boolean tumor;
    private final int chromosomeLength;

    public RefContextConsumer(boolean tumor, final int minQuality, @NotNull final GenomeRegion bounds,
            @NotNull final IndexedFastaSequenceFile refGenome, @NotNull final RefContextCandidates candidates) {
        this.bounds = bounds;
        this.refGenome = refGenome;
        this.minQuality = minQuality;
        this.candidates = candidates;
        this.tumor = tumor;
        chromosomeLength = refGenome.getSequence(bounds.chromosome()).length();

    }

    @Override
    public void accept(final SAMRecord record) {

        if (inBounds(record)) {

            int alignmentStart = record.getAlignmentStart();
            int alignmentEnd = record.getAlignmentEnd();

            if (record.getMappingQuality() >= minQuality) {

                byte[] refBases = refGenome.getSubsequenceAt(record.getContig(), alignmentStart, alignmentEnd).getBases();

                int readIndex = 0;
                int refBase = record.getAlignmentStart();

                for (final CigarElement e : record.getCigar().getCigarElements()) {
                    final int length = e.getLength();

                    switch (e.getOperator()) {
                        case M:
                        case EQ:
                        case X:
                            processAlignment(record, readIndex, refBase, length, refBases);
                            break;
                        case D:
                            //                            processDelete(record, length, readIndex - 1, refBase - 1, refBases);
                            break;
                    }

                    if (e.getOperator().consumesReferenceBases()) {
                        refBase += e.getLength();
                    }

                    if (e.getOperator().consumesReadBases()) {
                        readIndex += e.getLength();
                    }
                }

            } else {

                processSubprime(record);
            }

        }
    }

    private void processDelete(@NotNull final SAMRecord record, int length, int readIndex, int refPosition, byte[] refBases) {

        int refIndex = refPosition - record.getAlignmentStart();

        final String alt = new String(refBases, refIndex, 1);
        final String ref = new String(refBases, refIndex, length + 1);

        final RefContext refContext = candidates.refContext(record.getContig(), refPosition);
        if (refContext != null) {
            if (tumor) {
                refContext.altRead(ref, alt, new ReadContext(readIndex, record, refBases));
            } else {
                refContext.altRead(ref, alt);
            }
        }
    }

    private void processAlignment(@NotNull final SAMRecord record, int readBasesStartIndex, int refPositionStart, int alignmentLength,
            byte[] refBases) {

        int refBasesStartIndex = refPositionStart - record.getAlignmentStart();

        for (int i = 0; i < alignmentLength; i++) {

            long refPosition = refPositionStart + i;
            int readBaseIndex = readBasesStartIndex + i;
            int refBaseIndex = refBasesStartIndex + i;

            if (!inBounds(refPosition)) {
                continue;
            }

            final byte refByte = refBases[refBaseIndex];
            final String ref = String.valueOf((char) refByte);
            final byte readByte = record.getReadBases()[readBaseIndex];
            final int baseQuality = record.getBaseQualities()[readBaseIndex];

            final RefContext refContext = candidates.refContext(record.getContig(), refPosition);
            if (refContext != null) {

                if (readByte != refByte) {
                    final String alt = String.valueOf((char) readByte);
                    if (tumor) {
                        refContext.altRead(ref, alt, new ReadContext(readBaseIndex, record, refBases));
                    } else {
                        refContext.altRead(ref, alt);
                    }

                } else {
                    refContext.refRead(record.getMappingQuality(), baseQuality);
                }
            }
        }

    }

    private void processSubprime(@NotNull final SAMRecord record) {
        int refStart = record.getAlignmentStart();
        int refEnd = record.getAlignmentEnd();

        for (int refPosition = refStart; refPosition <= refEnd; refPosition++) {

            final RefContext refContext = candidates.refContext(record.getContig(), refPosition);
            if (refContext != null) {
                refContext.subprimeRead(record.getMappingQuality());
            }
        }
    }

    private boolean inBounds(final SAMRecord record) {
        return record.getEnd() >= bounds.start() && record.getStart() <= bounds.end();
    }

    private boolean inBounds(final long position) {
        return position >= bounds.start() && position <= bounds.end();
    }

}
