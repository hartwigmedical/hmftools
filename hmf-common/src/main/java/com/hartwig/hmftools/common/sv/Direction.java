package com.hartwig.hmftools.common.sv;

public enum Direction
{
    FORWARDS(1),
    REVERSE(-1),
    UNKNOWN(0),
    ;

    public final int Step;

    Direction(final int step) {
        Step = step;
    }

    public Direction opposite()
    {
        switch(this)
        {
            case FORWARDS:
                return REVERSE;
            case REVERSE:
                return FORWARDS;
            default:
                return UNKNOWN;
        }
    }

    public static Direction parse(final String input) {
        switch(input) {
            case "FORWARDS":
            case "1":
                return Direction.FORWARDS;
            case "REVERSE":
            case "-1":
                return Direction.REVERSE;
            case "UNKNOWN":
            case "0":
                return Direction.UNKNOWN;
            default:
                throw new IllegalArgumentException("Cannot parse " + input + " as direction");
        }
    }
}
