package com.hartwig.hmftools.sage.context;

import java.util.function.Consumer;

import com.hartwig.hmftools.common.region.GenomeRegion;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class RefContextConsumer implements Consumer<SAMRecord> {

    private static final Logger LOGGER = LogManager.getLogger(RefContextConsumer.class);


    /*
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
    */

    private static final String DISTANCE_FROM_REF_TAG = "NM";
    private final int minQuality;

    private final GenomeRegion bounds;
    private final IndexedFastaSequenceFile refGenome;
    private final RefContextCandidates candidates;
    private final boolean doReadContext;

    public RefContextConsumer(boolean readContext, final int minQuality, @NotNull final GenomeRegion bounds,
            @NotNull final IndexedFastaSequenceFile refGenome, @NotNull final RefContextCandidates candidates) {
        this.bounds = bounds;
        this.refGenome = refGenome;
        this.minQuality = minQuality;
        this.candidates = candidates;
        this.doReadContext = readContext;
    }

    @Override
    public void accept(final SAMRecord record) {

        if (inBounds(record)) {

            if (record.getMappingQuality() >= minQuality) {
                record.getAlignmentBlocks().forEach(x -> processPrimeAlignment(record, x));
            } else {
                record.getAlignmentBlocks().forEach(x -> processSubprimeAlignment(record, x));
            }

        }
    }

    private void processSubprimeAlignment(@NotNull final SAMRecord record, @NotNull final AlignmentBlock alignmentBlock) {
        long refStart = alignmentBlock.getReferenceStart();
        long refEnd = refStart + alignmentBlock.getLength() - 1;
        byte[] refBytes = refGenome.getSubsequenceAt(record.getContig(), refStart, refEnd).getBases();

        for (int refBytePosition = 0; refBytePosition < alignmentBlock.getLength(); refBytePosition++) {
            long position = refStart + refBytePosition;

            final byte refByte = refBytes[refBytePosition];
            final String ref = String.valueOf((char) refByte);

            RefContext refContext = candidates.refContext(record.getContig(), position, ref);
            if (refContext != null) {
                refContext.subprimeRead(record.getMappingQuality());
            }
        }
    }

    private void processPrimeAlignment(@NotNull final SAMRecord record, @NotNull final AlignmentBlock alignmentBlock) {

        long refStart = alignmentBlock.getReferenceStart();
        long refEnd = refStart + alignmentBlock.getLength() - 1;
        byte[] refBytes = refGenome.getSubsequenceAt(record.getContig(), refStart, refEnd).getBases();

        int readStart = alignmentBlock.getReadStart() - 1;

        for (int refBytePosition = 0; refBytePosition < alignmentBlock.getLength(); refBytePosition++) {

            long position = refStart + refBytePosition;
            int readBytePosition = readStart + refBytePosition;

            if (!inBounds(position)) {
                continue;
            }

            final byte refByte = refBytes[refBytePosition];
            final String ref = String.valueOf((char) refByte);
            final byte readByte = record.getReadBases()[readBytePosition];
            final int baseQuality = record.getBaseQualities()[readBytePosition];

            final RefContext refContext = candidates.refContext(record.getContig(), position, ref);
            if (refContext != null) {

                if (readByte != refByte) {
                    final String alt = String.valueOf((char) readByte);

                    long distanceFromAlignmentStart = position - record.getAlignmentStart();
                    long distanceFromAlignmentEnd = record.getAlignmentEnd() - position;
                    int minDistanceFromAlignment = (int) Math.min(distanceFromAlignmentStart, distanceFromAlignmentEnd);

                    if (doReadContext) {
                        refContext.altRead(alt,
                                record.getMappingQuality(),
                                baseQuality,
                                readBytePosition,
                                minDistanceFromAlignment,
                                new ReadContext(readBytePosition, record.getReadBases()));
                    } else {
                        refContext.altRead(alt,
                                record.getMappingQuality(),
                                baseQuality,
                                readBytePosition,
                                minDistanceFromAlignment);
                    }

                } else {
                    refContext.refRead(record.getMappingQuality(), baseQuality);
                }
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
