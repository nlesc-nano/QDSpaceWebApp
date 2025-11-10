# Single-stage on Lambda base (Amazon Linux 2)
FROM public.ecr.aws/lambda/python:3.11

# Set working directory
WORKDIR ${LAMBDA_TASK_ROOT}

# Install system deps with yum (for git, compilation, RDKit/cairo/pango deps)
RUN yum update -y && \
    yum install -y git gcc gcc-c++ make boost-devel cairo-devel pango-devel pkgconfig-devel libffi-devel openssl-devel mesa-libGL-devel && \
    yum clean all && \
    rm -rf /var/cache/yum

# Copy requirements.txt first (for caching)
COPY requirements.txt ${LAMBDA_TASK_ROOT}

# Install Python deps (includes Mangum, python-multipart, Git repos)
RUN pip install --no-cache-dir -r requirements.txt

# Copy the entire repo (including /docs for static, backend/, main.py)
COPY . ${LAMBDA_TASK_ROOT}

# Ensure /tmp for temp files
ENV TMPDIR=/tmp

# Lambda entry: Mangum handler
CMD ["main.handler"]