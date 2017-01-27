package com.hartwig.hmftools.retentionchecker;

class FastqHeaderKey {

    private final int tileNumber;
    private final int xCoordinate;
    private final int yCoordinate;

    FastqHeaderKey(final int tileNumber, final int xCoordinate, final int yCoordinate) {
        this.tileNumber = tileNumber;
        this.xCoordinate = xCoordinate;
        this.yCoordinate = yCoordinate;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        FastqHeaderKey that = (FastqHeaderKey) o;

        if (tileNumber != that.tileNumber) return false;
        if (xCoordinate == that.xCoordinate) if (yCoordinate == that.yCoordinate) return true;
        return false;
    }

    @Override
    public int hashCode() {
        int result = tileNumber;
        result = 31 * result + xCoordinate;
        result = 31 * result + yCoordinate;
        return result;
    }
}
