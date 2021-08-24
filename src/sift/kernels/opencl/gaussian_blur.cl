

__kernel void gaussian_blur(
    __global unsigned char *pixels,
    __global unsigned char *out,
    __global double* gaussian,
    const unsigned long heigth,
    const unsigned long width,
    const unsigned long window_width,
    const unsigned long chans
)
{
	int i, j, x, y, k, rows, cols, idx;
	
	float acc = 0;

	idx  = get_global_id(0);

	rows = (int)heigth;
	cols = (int)width;

	k    = window_width;
	y    = idx / cols;
	x    = idx % cols;

	for (int t = 0; t < chans - 1; t++) {

		int colStart = max(x - k/2, 0);
		int colEnd   = min(x + k/2, cols);

		int rowStart = max(y - k/2, 0);
		int rowEnd   = min(y + k/2, rows);

		int marginLeft   = clamp(x - k/2, -k, 0);
		int marginTop    = clamp(y - k/2, -k, 0);

		for (int r = rowStart; r < rowEnd; r++)
		{
			for (int c = colStart; c < colEnd; c++)
			{
				int gIdx = (r - rowStart - marginTop) * k + (c - colStart - marginLeft);
				int pIdx = (r * width + c) * chans + t;

				acc += gaussian[gIdx] * pixels[pIdx];
			}
		}
		
		out[(y * width + x) * chans + t] = (unsigned char)acc;
	}

	out[(y * cols + x) * chans + chans - 1] = 255;
}
