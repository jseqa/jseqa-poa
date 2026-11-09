package jseqa.poa;

import static org.junit.Assert.assertEquals;
import org.junit.Test;
import jseqa.poa.POGPtr.PointerType;

public class AlignPOGPtrMatrixPrimitiveTest {

    @Test
    public void uninitializedCellsRetainLegacyUndefinedPointer() {
        POGPtr pointer = new AlignPOGPtrMatrix(2, 2).get(1, 1);
        assertEquals(-1, pointer.ptr().first().intValue());
        assertEquals(-1, pointer.ptr().second().intValue());
        assertEquals(0, pointer.score());
        assertEquals(PointerType.UNDEF, pointer.type());
    }

    @Test
    public void storesPrimitivePointerComponents() {
        AlignPOGPtrMatrix matrix = new AlignPOGPtrMatrix(3, 4);
        matrix.set(2, 3, 1, 2, 17, PointerType.MIS);
        assertEquals(1, matrix.previousRow(2, 3));
        assertEquals(2, matrix.previousColumn(2, 3));
        assertEquals(17, matrix.score(2, 3));
        assertEquals(PointerType.MIS, matrix.type(2, 3));
    }

    @Test
    public void compatibilityGetterCreatesEquivalentPointer() {
        AlignPOGPtrMatrix matrix = new AlignPOGPtrMatrix(2, 2);
        matrix.set(1, 1, 0, 0, 5, PointerType.MATCH);
        POGPtr pointer = matrix.get(1, 1);
        assertEquals(0, pointer.ptr().first().intValue());
        assertEquals(0, pointer.ptr().second().intValue());
        assertEquals(5, pointer.score());
        assertEquals(PointerType.MATCH, pointer.type());
    }

    @Test
    public void resizePreservesObjectAndUpdatesActiveDimensions() {
        AlignPOGPtrMatrix matrix = new AlignPOGPtrMatrix(2, 8);
        matrix.resize(2, 3);
        assertEquals(2, matrix.nrows());
        assertEquals(3, matrix.ncols());
        matrix.resize(5, 12);
        assertEquals(5, matrix.nrows());
        assertEquals(12, matrix.ncols());
    }
}
