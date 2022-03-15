package com.hartwig.hmftools.serve.extraction.hotspot.tools;

import java.io.IOException;
import java.util.List;
import java.util.Set;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.codon.AminoAcids;
import com.hartwig.hmftools.common.variant.impact.VariantTranscriptImpact;
import com.hartwig.hmftools.serve.extraction.util.VCFWriterFactory;

import org.apache.commons.compress.utils.Lists;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.config.Configurator;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;

public class AnnotatedHotspotVCFCheckerPAVE {

    private static final Logger LOGGER = LogManager.getLogger(AnnotatedHotspotVCFCheckerPAVE.class);
    private static final boolean LOG_DEBUG = false;

    // Should be enabled when we know the hotspots have been lifted-over
    private static final boolean ENABLE_TRANSCRIPT_REF_GENOME_CURATION = false;

    private static final String NO_INPUT_PROTEIN = "-";

    private final Set<String> curatedTranscripts = Sets.newHashSet();

    public static void main(String[] args) throws IOException {
        if (LOG_DEBUG) {
            Configurator.setRootLevel(Level.DEBUG);
        }

        String annotatedHotspotVcf = System.getProperty("user.home") + "/hmf/tmp/KnownHotspots.somatic.37.pave.vcf";
        new AnnotatedHotspotVCFCheckerPAVE().run(annotatedHotspotVcf);
    }

    public void run(@NotNull String annotatedVcfFilePath) throws IOException {
        int totalCount = 0;
        int matchCount = 0;
        int approximateMatchCount = 0;
        int transcriptChangeLiftoverCount = 0;
        int diffCount = 0;

        LOGGER.info("Loading hotspots from '{}'", annotatedVcfFilePath);
        AbstractFeatureReader<VariantContext, LineIterator> reader =
                AbstractFeatureReader.getFeatureReader(annotatedVcfFilePath, new VCFCodec(), false);
        for (VariantContext variant : reader.iterator()) {
            totalCount++;

            String[] inputParts = variant.getAttributeAsString(VCFWriterFactory.INPUT_FIELD, Strings.EMPTY).split("\\|");
            String inputGene = inputParts[0];
            String inputTranscript = inputParts[1].equals("null") ? null : inputParts[1];
            String inputProteinAnnotation = inputParts[2];

            String formattedHotspot = formatHotspot(variant);
            if (inputProteinAnnotation.equals(NO_INPUT_PROTEIN)) {
                LOGGER.debug("Skipping non-coding hotspot on '{}'", formattedHotspot);
                matchCount++;
            } else {
                List<VariantTranscriptImpact> annotations = toVariantTranscriptImpact(variant.getAttributeAsStringList("PAVE_TI", ""));

                MatchType match = determineMatch(inputTranscript, inputProteinAnnotation, annotations);

                switch (match) {
                    case IDENTICAL: {
                        LOGGER.debug("Found a match amongst candidate transcripts for '{}' on '{}", inputProteinAnnotation, inputGene);
                        matchCount++;
                        break;
                    }
                    case APPROXIMATE: {
                        LOGGER.debug("Found approximate match amongst candidate transcripts for '{}' on '{}",
                                inputProteinAnnotation,
                                inputGene);
                        matchCount++;
                        approximateMatchCount++;
                        break;
                    }
                    case LIFTOVER_RETIRED_OR_CHANGED_TRANSCRIPT: {
                        LOGGER.debug("Match considered valid because of retired or changed transcript during liftover");
                        transcriptChangeLiftoverCount++;
                        matchCount++;
                        break;
                    }
                    case NO_MATCH: {
                        //                        LOGGER.warn("Could not match inputTranscript {}, input protein {} for input gene {}",
                        //                                inputTranscript,
                        //                                inputProteinAnnotation,
                        //                                inputGene);

                        for (VariantTranscriptImpact variantTranscriptImpact : annotations) {
                            if (variantTranscriptImpact.Transcript.equals(inputTranscript)) {

                                LOGGER.warn(
                                        "Could not match inputTranscript {}, paveTranscript {}, input protein {}, pave protein {}, pave coding {} for input gene {}, pave gene {}",
                                        inputTranscript,
                                        variantTranscriptImpact.Transcript,
                                        inputProteinAnnotation,
                                        variantTranscriptImpact.HgvsProtein,
                                        variantTranscriptImpact.HgvsCoding,
                                        inputGene,
                                        variantTranscriptImpact.GeneName);

                            }
                        }

                        diffCount++;

                        break;
                    }
                }
            }
        }

        LOGGER.info("Done comparing {} records: {} matches (of which {} are approximate and {} are due to transcript liftover changes)"
                + " and {} differences found.", totalCount, matchCount, approximateMatchCount, transcriptChangeLiftoverCount, diffCount);

        // current result will be
        // Done comparing 9452 records: 9452 matches (of which 30 are approximate and 0 are due to transcript liftover changes) and 0 differences found.
        checkForUnusedMappings();
    }

    @NotNull
    private static List<VariantTranscriptImpact> toVariantTranscriptImpact(@NotNull final List<String> annotation) {
        List<VariantTranscriptImpact> variantTranscriptImpactList = Lists.newArrayList();
        for (int i = 0; i <= annotation.size() - 1; i++) {
            variantTranscriptImpactList.add(VariantTranscriptImpact.fromVcfData(annotation.get(i)));
        }

        return variantTranscriptImpactList;
    }

    @NotNull
    private static String formatHotspot(@NotNull VariantContext variant) {
        return variant.getContig() + ":" + variant.getStart() + " " + variant.getReference().getBaseString() + ">"
                + variant.getAlternateAllele(0).getBaseString();
    }

    @NotNull
    private MatchType determineMatch(@Nullable String inputTranscript, @NotNull String inputProteinAnnotation,
            @NotNull List<VariantTranscriptImpact> annotations) {
        VariantTranscriptImpact specificAnnotation = annotationForTranscript(annotations, inputTranscript);
        if (specificAnnotation != null) {
            return matchOnSpecificAnnotation(inputTranscript, inputProteinAnnotation, specificAnnotation);
        } else {
            // In case input transcript is missing or can't be found, we try to match against any transcript.
            // This could be tricky in case a variant was generated from 37 and is now being evaluated on 38 with different transcript IDs.
            return matchOnAnyTranscript(inputProteinAnnotation, annotations);
        }
    }

    @NotNull
    private MatchType matchOnSpecificAnnotation(@NotNull String inputTranscript, @NotNull String inputProteinAnnotation,
            @NotNull VariantTranscriptImpact specificAnnotation) {
        String PaveProteinAnnotation = AminoAcids.forceSingleLetterProteinAnnotation(specificAnnotation.HgvsProtein);
        return matchAnnotation(inputTranscript, inputProteinAnnotation, PaveProteinAnnotation, specificAnnotation.Effects);
    }

    @NotNull
    private MatchType matchOnAnyTranscript(@NotNull String inputProteinAnnotation, @NotNull List<VariantTranscriptImpact> annotations) {
        MatchType matchedMatchType = MatchType.NO_MATCH;
        for (VariantTranscriptImpact annotation : annotations) {
            // We only want to consider transcript features with coding impact.
            if (!annotation.HgvsProtein.isEmpty()) {
                String PaveProteinAnnotation = AminoAcids.forceSingleLetterProteinAnnotation(annotation.HgvsProtein);
                MatchType match = matchAnnotation(annotation.Transcript, inputProteinAnnotation, PaveProteinAnnotation, annotation.Effects);
                if (match != MatchType.NO_MATCH && matchedMatchType == MatchType.NO_MATCH) {
                    matchedMatchType = match;
                }
            }
        }

        return matchedMatchType;
    }

    @Nullable
    private static VariantTranscriptImpact annotationForTranscript(@NotNull List<VariantTranscriptImpact> annotations,
            @Nullable String transcript) {
        for (VariantTranscriptImpact annotation : annotations) {
            if (annotation.Transcript.equals(transcript)) {
                return annotation;
            }
        }
        return null;
    }

    @NotNull
    private MatchType matchAnnotation(@NotNull String transcript, @NotNull String inputAnnotation, @NotNull String PaveProteinAnnotation,
            @NotNull String effect) {
        String curatedInputAnnotation = curateStartCodonAnnotation(inputAnnotation);
        curatedInputAnnotation = curateSynonymous(curatedInputAnnotation, effect, PaveProteinAnnotation);
        if (curatedInputAnnotation.equals(PaveProteinAnnotation)) {
            return MatchType.IDENTICAL;
        }

        if (isApproximateIndelMatch(inputAnnotation, PaveProteinAnnotation)) {
            return MatchType.APPROXIMATE;
        }

        if (ENABLE_TRANSCRIPT_REF_GENOME_CURATION && (retiredTranscriptCheck(transcript, PaveProteinAnnotation) || changedTranscriptCheck(
                transcript,
                inputAnnotation,
                PaveProteinAnnotation))) {
            return MatchType.LIFTOVER_RETIRED_OR_CHANGED_TRANSCRIPT;
        }

        return MatchType.NO_MATCH;
    }

    @NotNull
    private static String curateStartCodonAnnotation(@NotNull String serveAnnotation) {
        if (serveAnnotation.startsWith("p.M1") && serveAnnotation.length() == 5) {
            return "p.M1?";
        } else {
            return serveAnnotation;
        }
    }

    @NotNull
    private static String curateSynonymous(@NotNull String serveAnnotation, @NotNull String effect, @NotNull String PaveProteinAnnotation) {
        if (effect.contains("synonymous_variant") || effect.contains("start_lost") || effect.contains("5_prime_UTR_variant")) {
            return PaveProteinAnnotation;
        } else {
            return serveAnnotation;
        }
    }

    @VisibleForTesting
    static boolean isApproximateIndelMatch(@NotNull String inputAnnotation, @NotNull String snpEffAnnotation) {
        if (inputAnnotation.contains("del") || inputAnnotation.contains("ins") || inputAnnotation.contains("dup")) {
            int inputStartCodon = extractDigitWithIndex(inputAnnotation, 1);
            Integer inputEndCodon = extractDigitWithIndex(inputAnnotation, 2);
            int snpeffStartCodon = extractDigitWithIndex(snpEffAnnotation, 1);

            int maxDistance = inputEndCodon != null ? 1 + inputEndCodon - inputStartCodon : 3;

            boolean aminoAcidCheck = true;
            if (inputAnnotation.endsWith("del")) {
                // For deletes we allow for difference in alignment since they can often be realigned.
                maxDistance = 20;
                if (!inputAnnotation.contains("_")) {
                    // Single AA deletes have to match on specific AA
                    aminoAcidCheck = inputAnnotation.substring(2, 3).equals(snpEffAnnotation.substring(2, 3));
                }
            }

            return aminoAcidCheck && Math.abs(inputStartCodon - snpeffStartCodon) <= maxDistance;
        }

        return false;
    }

    @Nullable
    private static Integer extractDigitWithIndex(@NotNull String string, int digitIndex) {
        StringBuilder digitBuilder = new StringBuilder();
        boolean isExtracting = false;
        int numExtracted = 0;
        for (int i = 0; i < string.length(); i++) {
            char letter = string.charAt(i);
            if (isInteger(letter)) {
                if (!isExtracting) {
                    isExtracting = true;
                    numExtracted++;
                }

                if (isExtracting && numExtracted == digitIndex) {
                    digitBuilder.append(letter);
                }
            } else {
                isExtracting = false;
            }
        }

        String digitString = digitBuilder.toString();
        return !digitString.isEmpty() ? Integer.parseInt(digitString) : null;
    }

    private static boolean isInteger(char letter) {
        try {
            Integer.parseInt(String.valueOf(letter));
            return true;
        } catch (NumberFormatException exception) {
            return false;
        }
    }

    private boolean retiredTranscriptCheck(@NotNull String transcript, @NotNull String snpeffAnnotation) {
        if (AnnotatedHotspotCurationFactory.RETIRED_TRANSCRIPTS.contains(transcript) && snpeffAnnotation.isEmpty()) {
            // In case we know a transcript has been retired from coding duty in certain ref genomes we accept the diff when empty.
            curatedTranscripts.add(transcript);
            return true;
        }

        return false;
    }

    private boolean changedTranscriptCheck(@NotNull String transcript, @NotNull String inputAnnotation, @NotNull String snpeffAnnotation) {
        if (AnnotatedHotspotCurationFactory.CHANGED_TRANSCRIPTS.contains(transcript)) {
            // In case transcripts have different versions across ref genomes we assume the AA change is the same, just a different position
            // Eg p.PxxS should match for any xx
            curatedTranscripts.add(transcript);
            return inputAnnotation.substring(0, 3).equals(snpeffAnnotation.substring(0, 3))
                    && inputAnnotation.substring(inputAnnotation.length()).equals(snpeffAnnotation.substring(snpeffAnnotation.length()));
        }

        return false;
    }

    private void checkForUnusedMappings() {
        if (ENABLE_TRANSCRIPT_REF_GENOME_CURATION) {
            int unusedTranscriptCurationCount = 0;
            for (String transcript : AnnotatedHotspotCurationFactory.RETIRED_TRANSCRIPTS) {
                if (!curatedTranscripts.contains(transcript)) {
                    LOGGER.warn("Transcript '{}' labeled as 'retired' has not been used in curation", transcript);
                    unusedTranscriptCurationCount++;
                }
            }

            for (String transcript : AnnotatedHotspotCurationFactory.CHANGED_TRANSCRIPTS) {
                if (!curatedTranscripts.contains(transcript)) {
                    LOGGER.warn("Transcript '{}' labeled as 'changed' has not been used in curation", transcript);
                    unusedTranscriptCurationCount++;
                }
            }

            LOGGER.info("Analyzed usage of transcript curation. Found {} unused curation keys.", unusedTranscriptCurationCount);
        }
    }

    private enum MatchType {
        IDENTICAL,
        APPROXIMATE,
        LIFTOVER_RETIRED_OR_CHANGED_TRANSCRIPT,
        NO_MATCH
    }
}
