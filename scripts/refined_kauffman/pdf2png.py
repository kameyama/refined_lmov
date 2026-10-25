"""CoreGraphics を ctypes で呼んで PDF のページを PNG にする（poppler 無しの環境用）。使い方：pdf2png.py in.pdf outdir first last scale"""
import ctypes, ctypes.util, sys, zlib, struct, os
cg = ctypes.CDLL("/System/Library/Frameworks/CoreGraphics.framework/CoreGraphics")
cf = ctypes.CDLL("/System/Library/Frameworks/CoreFoundation.framework/CoreFoundation")
class CGPoint(ctypes.Structure): _fields_ = [("x", ctypes.c_double), ("y", ctypes.c_double)]
class CGSize(ctypes.Structure): _fields_ = [("w", ctypes.c_double), ("h", ctypes.c_double)]
class CGRect(ctypes.Structure): _fields_ = [("origin", CGPoint), ("size", CGSize)]
cf.CFURLCreateFromFileSystemRepresentation.restype = ctypes.c_void_p
cf.CFURLCreateFromFileSystemRepresentation.argtypes = [ctypes.c_void_p, ctypes.c_char_p, ctypes.c_long, ctypes.c_bool]
cg.CGPDFDocumentCreateWithURL.restype = ctypes.c_void_p; cg.CGPDFDocumentCreateWithURL.argtypes = [ctypes.c_void_p]
cg.CGPDFDocumentGetNumberOfPages.restype = ctypes.c_size_t; cg.CGPDFDocumentGetNumberOfPages.argtypes = [ctypes.c_void_p]
cg.CGPDFDocumentGetPage.restype = ctypes.c_void_p; cg.CGPDFDocumentGetPage.argtypes = [ctypes.c_void_p, ctypes.c_size_t]
cg.CGPDFPageGetBoxRect.restype = CGRect; cg.CGPDFPageGetBoxRect.argtypes = [ctypes.c_void_p, ctypes.c_int]
cg.CGColorSpaceCreateDeviceRGB.restype = ctypes.c_void_p
cg.CGBitmapContextCreate.restype = ctypes.c_void_p
cg.CGBitmapContextCreate.argtypes = [ctypes.c_void_p, ctypes.c_size_t, ctypes.c_size_t, ctypes.c_size_t, ctypes.c_size_t, ctypes.c_void_p, ctypes.c_uint32]
cg.CGContextSetRGBFillColor.argtypes = [ctypes.c_void_p, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double]
cg.CGContextFillRect.argtypes = [ctypes.c_void_p, CGRect]
cg.CGContextScaleCTM.argtypes = [ctypes.c_void_p, ctypes.c_double, ctypes.c_double]
cg.CGContextTranslateCTM.argtypes = [ctypes.c_void_p, ctypes.c_double, ctypes.c_double]
cg.CGContextDrawPDFPage.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
cg.CGBitmapContextGetData.restype = ctypes.c_void_p; cg.CGBitmapContextGetData.argtypes = [ctypes.c_void_p]
cg.CGContextRelease.argtypes = [ctypes.c_void_p]

def write_png(path, w, h, rgba):
    raw = bytearray()
    for y in range(h):
        row = rgba[y*w*4:(y+1)*w*4]
        raw.append(0)
        raw.extend(bytes(row[i] for i in range(len(row)) if i % 4 != 3))  # RGB
    def chunk(tag, data):
        c = struct.pack(">I", len(data)) + tag + data
        return c + struct.pack(">I", zlib.crc32(tag + data) & 0xffffffff)
    png = b"\x89PNG\r\n\x1a\n" + chunk(b"IHDR", struct.pack(">IIBBBBB", w, h, 8, 2, 0, 0, 0)) + chunk(b"IDAT", zlib.compress(bytes(raw), 6)) + chunk(b"IEND", b"")
    open(path, "wb").write(png)

path, outdir, first, last, scale = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), float(sys.argv[5])
url = cf.CFURLCreateFromFileSystemRepresentation(None, path.encode(), len(path.encode()), False)
doc = cg.CGPDFDocumentCreateWithURL(url)
n = cg.CGPDFDocumentGetNumberOfPages(doc)
print("pages:", n)
cs = cg.CGColorSpaceCreateDeviceRGB()
for p in range(first, min(last, n)+1):
    page = cg.CGPDFDocumentGetPage(doc, p)
    box = cg.CGPDFPageGetBoxRect(page, 0)
    w, h = int(box.size.w*scale), int(box.size.h*scale)
    ctx = cg.CGBitmapContextCreate(None, w, h, 8, w*4, cs, 5)  # kCGImageAlphaNoneSkipLast
    cg.CGContextSetRGBFillColor(ctx, 1, 1, 1, 1); cg.CGContextFillRect(ctx, CGRect(CGPoint(0, 0), CGSize(w, h)))
    cg.CGContextScaleCTM(ctx, scale, scale); cg.CGContextTranslateCTM(ctx, -box.origin.x, -box.origin.y)
    cg.CGContextDrawPDFPage(ctx, page)
    data = cg.CGBitmapContextGetData(ctx)
    buf = ctypes.string_at(data, w*h*4)
    # CoreGraphics の原点は左下なので行を反転
    rows = [buf[y*w*4:(y+1)*w*4] for y in range(h)]
    write_png(os.path.join(outdir, f"page{p:03d}.png"), w, h, b"".join(rows))
    cg.CGContextRelease(ctx)
    print("wrote page", p, w, h)
