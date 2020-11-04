package com.hartwig.hmftools.vicc.annotation;

import java.util.Map;
import java.util.Set;

import com.google.common.collect.Maps;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class ComplexClassifier {

    private static final Map<String, Set<String>> COMPLEX_EVENTS_PER_GENE = Maps.newHashMap();

    static {
        //        Set<String> alkSet = Sets.newHashSet("ALK inframe insertion (1151T)");
        //        COMPLEX_EVENTS_PER_GENE.put("ALK", alkSet);
        //
        //        Set<String> apcSet = Sets.newHashSet("APC p.I1557*fs*1", "APC p.K1310*fs*1", "APC p.G1466*fs*1");
        //        COMPLEX_EVENTS_PER_GENE.put("APC", apcSet);
        //
        //        Set<String> arSet = Sets.newHashSet("ARv567es", "SPLICE VARIANT 7");
        //        COMPLEX_EVENTS_PER_GENE.put("AR", arSet);
        //
        //        Set<String> b2mSet = Sets.newHashSet("B2M .");
        //        COMPLEX_EVENTS_PER_GENE.put("B2M", b2mSet);
        //
        //        Set<String> brafSet = Sets.newHashSet("DEL 485-490", "L485_P490>Y", ); COMPLEX_EVENTS_PER_GENE.put("BRAF", set);
        //
        //        Set<String> brca1Set = Sets.newHashSet("BRCA1 L392Qfs*5", "BRCA1 .");
        //        COMPLEX_EVENTS_PER_GENE.put("BRCA1", set);
        //
        //        Set<String> brca2Set = Sets.newHashSet("");
        //        COMPLEX_EVENTS_PER_GENE.put("BRCA2", set);
        //
        //        Set<String> btkSet = Sets.newHashSet("");
        //        COMPLEX_EVENTS_PER_GENE.put("BTK", set);
        //
        //        Set<String> ccnd1Set = Sets.newHashSet("");
        //        COMPLEX_EVENTS_PER_GENE.put("CCND1", set);
        //
        //        Set<String> ccnd3Set = Sets.newHashSet("");
        //        COMPLEX_EVENTS_PER_GENE.put("CCND3", set);
        //
        //        Set<String> chek2Set = Sets.newHashSet("");
        //        COMPLEX_EVENTS_PER_GENE.put("CHEK2", set);
        //
        //        Set<String> Set = Sets.newHashSet("");
        //        COMPLEX_EVENTS_PER_GENE.put("", Set);
        //
        //        Set<String> Set = Sets.newHashSet("");
        //        COMPLEX_EVENTS_PER_GENE.put("", Set);

    }

    public static boolean isComplexEvent(@NotNull String featureName, @Nullable String gene) {
        Set<String> entriesForGene = COMPLEX_EVENTS_PER_GENE.get(gene);
        if (entriesForGene != null) {
            return entriesForGene.contains(featureName);
        }

        if (featureName.trim().endsWith(".")) {
            // Not sure what a dot means!
            return true;
        } else if (featureName.split("\\*").length > 2) {
            // Some hotspots contain multiple stop codons.
            return true;
        } else {
            String lowerFeature = featureName.toLowerCase();
            int fsLocation = lowerFeature.indexOf("fs");
            return fsLocation > 1 && !isInteger(lowerFeature.substring(fsLocation - 1, fsLocation));
        }
    }

    private static boolean isInteger(@NotNull String string) {
        try {
            Integer.parseInt(string);
            return true;
        } catch (NumberFormatException e) {
            return false;
        }
    }
}
