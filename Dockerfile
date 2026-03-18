# Lambda deployment (Amazon Linux 2023 base)
FROM public.ecr.aws/lambda/python:3.12

# Set working directory
WORKDIR ${LAMBDA_TASK_ROOT}

# Install system dependencies
RUN dnf update -y && \
    dnf install -y \
    gcc gcc-c++ make \
    python3-devel \
    boost-devel \
    cairo-devel pango-devel \
    pkgconfig \
    libffi-devel \
    mesa-libGL-devel && \
    dnf clean all && \
    rm -rf /var/cache/dnf

# Copy requirements.txt first (for caching)
COPY requirements.txt ${LAMBDA_TASK_ROOT}

# Install Python deps (includes Mangum, python-multipart, Git repos)
RUN pip install --no-cache-dir -r requirements.txt

# Copy application code
COPY . ${LAMBDA_TASK_ROOT}

# Set temp directory
ENV TMPDIR=/tmp

# Lambda handler
CMD ["main.handler"]