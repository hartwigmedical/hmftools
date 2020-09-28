package com.hartwig.hmftools.vicc.annotation;

import java.util.List;

import com.hartwig.hmftools.vicc.datamodel.ViccEntry;
import com.hartwig.hmftools.vicc.datamodel.cgi.Cgi;
import com.hartwig.hmftools.vicc.datamodel.civic.Civic;
import com.hartwig.hmftools.vicc.datamodel.oncokb.OncoKb;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class TranscriptExtractor {

    private TranscriptExtractor() {
    }

    @Nullable
    public static String extractTranscriptId(@NotNull ViccEntry viccEntry) {
        switch (viccEntry.source()) {
            case ONCOKB: {
                String transcriptId = ((OncoKb) viccEntry.kbSpecificObject()).transcriptID();
                return !transcriptId.isEmpty() ? transcriptId : null;
            }
            case CIVIC: {
                String transcriptId = ((Civic) viccEntry.kbSpecificObject()).coordinates().representativeTranscript();
                if (transcriptId == null || transcriptId.isEmpty()) {
                    return null;
                }
                int versionSeparator = transcriptId.indexOf(".");
                String transcript = versionSeparator >= 0 ? transcriptId.substring(0, versionSeparator) : transcriptId;
                // CIVIC occasionally just contains "ENST" or "ENST0000" as transcript. Need to filter that out.
                return transcript.length() == 15 ? transcript : null;
            }
            case CGI: {
                List<String> transcripts = ((Cgi) viccEntry.kbSpecificObject()).transcripts();
                // Even though a vicc entry could have multiple different transcripts, in practice every vicc entry with more than one
                // transcript only has one distinct transcript.
                String transcript = !transcripts.isEmpty() ? transcripts.get(0) : Strings.EMPTY;
                return !transcript.isEmpty() ? transcript : null;
            }
            case JAX:
            default:
                return null;
        }
    }
}
