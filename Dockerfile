# Lambda deployment (Amazon Linux 2 base)
FROM public.ecr.aws/lambda/python:3.11

# Set working directory
WORKDIR ${LAMBDA_TASK_ROOT}

# Install system dependencies
RUN yum update -y && \
    yum install -y git gcc gcc-c++ make boost-devel cairo-devel pango-devel pkgconfig-devel libffi-devel openssl-devel mesa-libGL-devel && \
    yum clean all && \
    rm -rf /var/cache/yum

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