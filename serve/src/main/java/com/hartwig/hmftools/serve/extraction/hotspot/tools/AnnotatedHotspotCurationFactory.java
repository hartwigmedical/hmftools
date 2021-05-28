package com.hartwig.hmftools.serve.extraction.hotspot.tools;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;

import org.jetbrains.annotations.NotNull;

final class AnnotatedHotspotCurationFactory {

    public static final Set<String> RETIRED_TRANSCRIPTS = createRetiredTranscriptWhitelist();
    public static final Set<String> CHANGED_TRANSCRIPTS = createChangedTranscriptWhitelist();
    public static final Map<String, Map<String, List<String>>> SERVE_TO_SNPEFF_MAPPINGS_PER_TRANSCRIPT = createTranscriptMappings();
    public static final Map<String, Map<String, List<String>>> SERVE_TO_SNPEFF_MAPPINGS_PER_GENE = createGeneMappings();

    @NotNull
    private static Set<String> createRetiredTranscriptWhitelist() {
        Set<String> transcripts = Sets.newHashSet();

        // This transcript on RHOA is coding in v37 but non-coding in v38
        transcripts.add("ENST00000431929");

        // This transcript on CDKN2A has retired in v38
        transcripts.add("ENST00000361570");

        // This transcript on BCL2L12 has completely changed in 38 and can't match anymore.
        transcripts.add("ENST00000246784");

        return transcripts;
    }

    @NotNull
    private static Set<String> createChangedTranscriptWhitelist() {
        Set<String> transcripts = Sets.newHashSet();

        // This MYC transcript is v2 in 37 and v6 in 38 and they differ by 15 AA.
        transcripts.add("ENST00000377970");

        // This TCF7L2 transcript in v1 in 37 and v5 in 38 and they differ by 5 AA.
        transcripts.add("ENST00000543371");

        // This FGFR2 transcript is v6 in 37 and v10 in 38 and they differ by 36 AA (spread across transcript).
        transcripts.add("ENST00000351936");

        // This MED12 transcript is v6 in 37 and v10 in 38 and they differ with 53 AA (spread across transcript)
        transcripts.add("ENST00000333646");

        return transcripts;
    }

    @NotNull
    private static Map<String, Map<String, List<String>>> createTranscriptMappings() {
        Map<String, Map<String, List<String>>> serveToSnpEffTranscriptMappings = Maps.newHashMap();

        serveToSnpEffTranscriptMappings.put("ENST00000288135", createKITTranscriptMap());
        serveToSnpEffTranscriptMappings.put("ENST00000275493", createEGFRTranscriptMap());
        serveToSnpEffTranscriptMappings.put("ENST00000269571", createERBB2TranscriptMap());
        serveToSnpEffTranscriptMappings.put("ENST00000263967", createPIK3CATranscriptMap());
        serveToSnpEffTranscriptMappings.put("ENST00000374690", createARTranscriptMap());
        serveToSnpEffTranscriptMappings.put("ENST00000231790", createMLH1TranscriptMap());

        return serveToSnpEffTranscriptMappings;
    }

    @NotNull
    private static Map<String, List<String>> createKITTranscriptMap() {
        Map<String, List<String>> map = Maps.newHashMap();
        map.put("p.S501_A502insAY", Lists.newArrayList("p.A502_Y503insYA"));
        map.put("p.V560del", Lists.newArrayList("p.V559del"));
        return map;
    }

    @NotNull
    private static Map<String, List<String>> createEGFRTranscriptMap() {
        Map<String, List<String>> map = Maps.newHashMap();
        map.put("p.V769_D770insASV", Lists.newArrayList("p.A767_V769dup"));
        return map;
    }

    @NotNull
    private static Map<String, List<String>> createERBB2TranscriptMap() {
        Map<String, List<String>> map = Maps.newHashMap();
        map.put("p.A771_Y772insYVMA", Lists.newArrayList("p.Y772_V773insVMAY"));
        map.put("p.Y772_A775dup", Lists.newArrayList("p.Y772_V773insVMAY"));
        map.put("p.M774_A775insAYVM", Lists.newArrayList("p.A775_G776insYVMA"));
        map.put("p.G778_P780dup", Lists.newArrayList("p.G778_S779insSPG"));
        return map;
    }

    @NotNull
    private static Map<String, List<String>> createPIK3CATranscriptMap() {
        Map<String, List<String>> map = Maps.newHashMap();
        map.put("p.E109del", Lists.newArrayList("p.E110del"));
        return map;
    }

    @NotNull
    private static Map<String, List<String>> createARTranscriptMap() {
        Map<String, List<String>> map = Maps.newHashMap();
        map.put("p.Q86del", Lists.newArrayList("p.Q87del", "p.Q88del", "p.Q89del", "p.Q90del", "p.Q91del"));
        return map;
    }

    @NotNull
    private static Map<String, List<String>> createMLH1TranscriptMap() {
        Map<String, List<String>> map = Maps.newHashMap();
        map.put("p.K618del", Lists.newArrayList("p.K616del", "p.K617del"));
        return map;
    }

    @NotNull
    private static Map<String, Map<String, List<String>>> createGeneMappings() {
        Map<String, Map<String, List<String>>> serveToSnpEffGeneMappings = Maps.newHashMap();
        serveToSnpEffGeneMappings.put("CDK11A", createCDK11AGeneMap());
        return serveToSnpEffGeneMappings;
    }

    @NotNull
    private static Map<String, List<String>> createCDK11AGeneMap() {
        Map<String, List<String>> map = Maps.newHashMap();
        map.put("p.K115_E116dup", Lists.newArrayList("p.E116_R117insKE"));
        return map;
    }

    private AnnotatedHotspotCurationFactory() {
    }
}
