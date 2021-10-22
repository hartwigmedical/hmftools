package com.hartwig.hmftools.patientreporter.cfreport.components;

import com.hartwig.hmftools.patientreporter.cfreport.ReportResources;
import com.itextpdf.kernel.colors.DeviceRgb;
import com.itextpdf.layout.element.Text;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public class Icon {

    public enum IconType {
        LEVEL_A,
        LEVEL_B,
        LEVEL_C,
        LEVEL_D,
        MATCH_BROAD,
        MATCH_SPECIFIC,
        TREATMENT,
        INVALID;

        @NotNull
        String text() {
            switch (this) {
                case LEVEL_A:
                    return "A";
                case LEVEL_B:
                    return "B";
                case LEVEL_C:
                    return "C";
                case LEVEL_D:
                    return "D";
                case MATCH_BROAD:
                    return "\u00a4";
                case MATCH_SPECIFIC:
                    return "\u00b0";
                case TREATMENT:
                    return "\u00b1";
                default:
                    return "";
            }
        }

        DeviceRgb defaultColor() {
            switch (this) {
                case LEVEL_A:
                    return ReportResources.PALETTE_CYAN;
                case LEVEL_B:
                    return ReportResources.PALETTE_BLUE;
                case LEVEL_C:
                    return ReportResources.PALETTE_DARK_BLUE;
                case LEVEL_D:
                    return ReportResources.PALETTE_VIOLET;
                case MATCH_BROAD:
                case MATCH_SPECIFIC:
                    return ReportResources.PALETTE_MID_BLUE;
                default:
                    return ReportResources.PALETTE_LIGHT_GREY;
            }
        }
    }

    private static final DeviceRgb[] TREATMENT_PALETTE =
            { new DeviceRgb(110, 197, 213)};

    @NotNull
    public static Text createLevelIcon(@NotNull String level) {
        IconType iconType;
        switch (level.toUpperCase()) {
            case "A":
                iconType = IconType.LEVEL_A;
                break;
            case "B":
                iconType = IconType.LEVEL_B;
                break;
            case "C":
                iconType = IconType.LEVEL_C;
                break;
            case "D":
                iconType = IconType.LEVEL_D;
                break;
            default:
                return new Text(Strings.EMPTY);
        }

        return createIcon(iconType);
    }

    @NotNull
    public static Text createTreatmentIcon(@NotNull String treatmentName) {
        int charCode = !treatmentName.isEmpty() ? (int) treatmentName.charAt(0) : 0;
        int colorIndex = charCode % TREATMENT_PALETTE.length;

        return createIcon(IconType.TREATMENT, TREATMENT_PALETTE[colorIndex]);
    }

    @NotNull
    public static Text createIcon(@NotNull IconType iconType) {
        return createIcon(iconType, iconType.defaultColor());
    }

    @NotNull
    private static Text createIcon(@NotNull IconType iconType, @NotNull DeviceRgb color) {
        if (iconType == IconType.INVALID) {
            return new Text("");
        }

        return new Text(iconType.text()).setFont(ReportResources.iconFont()).setFontColor(color).setCharacterSpacing(1);
    }
}
